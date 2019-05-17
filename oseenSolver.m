function [Uwn,pOld,tOld] = OseenSolver(p,t,f,fb)


waveNumbers=-1:2;
maxWaveNum=length(waveNumbers);
b=boundedges(p,t);
Dirichlet_e=b;
p=p';

ubar =@(x) 1+0*norm(x)^2;
pOld = p;
tOld = t;

% Quadratic elements contain six nodes per triangle: add three nodes at the middle of the edges of each triangle
numberOfNodes.old = size(p, 2);
numberOfElements = size(t, 1);
S = zeros(numberOfNodes.old);
counter = numberOfNodes.old + 1;
 
for e = 1:numberOfElements
    nodes = t(e, :); 
    if (S(nodes(1), nodes(2)) == 0)
        S(nodes(1), nodes(2)) = counter;
        S(nodes(2), nodes(1)) = counter;
        p(:, counter) = mean(p(:, [nodes(1) nodes(2)]), 2);
        counter = counter + 1;
    end
    if (S(nodes(2), nodes(3)) == 0)
        S(nodes(2), nodes(3)) = counter;
        S(nodes(3), nodes(2)) = counter;
        p(:, counter) = mean(p(:, [nodes(2) nodes(3)]), 2);
        counter = counter + 1;
    end
    if (S(nodes(1), nodes(3)) == 0)
        S(nodes(1), nodes(3)) = counter;
        S(nodes(3), nodes(1)) = counter;
        p(:, counter) = mean(p(:, [nodes(1) nodes(3)]), 2);
        counter = counter + 1;
    end
    t(e, 4) = S(nodes(1), nodes(2));
    t(e, 5) = S(nodes(2), nodes(3));
    t(e, 6) = S(nodes(1), nodes(3));
end
numberOfNodes.new = size(p, 2);
wnl=3*numberOfNodes.new+numberOfNodes.old;
Dirichlet = [];
for i=1 : length(Dirichlet_e)
    nodes = Dirichlet_e(i, :);
    Dirichlet = [Dirichlet; [nodes(1), S(nodes(1), nodes(2)), nodes(2)]];
end
Dirichlet=unique(Dirichlet);
Fb=zeros(length(Dirichlet)*maxWaveNum,1);
for i=1:maxWaveNum
    waveNum=waveNumbers(i);
 Fb((1:3*length(Dirichlet))+(i-1)*3*length(Dirichlet))=reshape(fb(p(:,Dirichlet)',waveNum),3*length(Dirichlet),1);
end

% Initialisation of K, F
K = sparse(numberOfNodes.new,numberOfNodes.new);
M = sparse(numberOfNodes.new,numberOfNodes.new);
B = sparse(2*numberOfNodes.new,numberOfNodes.old);
B3 = sparse(numberOfNodes.new,numberOfNodes.old);
S1 = sparse(numberOfNodes.new, numberOfNodes.new);
S2 = sparse(numberOfNodes.new, numberOfNodes.new);
T = sparse(numberOfNodes.new, numberOfNodes.new);
F = zeros(wnl*maxWaveNum,1);
 
% Gaussian quadrature points & weights
eta_xi = [2/3 1/6 1/6;
          1/6 2/3 1/6;
          1/6 1/6 2/3];
wOmega = [1/3 1/3 1/3];
 
% Assembly of K and F
set(gcf, 'color', 'w');
hold on
for e = 1:numberOfElements
    nodes = t(e, :);
    
    % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
    P = [ones(1, 6); 
         p(:, nodes); 
         p(1, nodes).^2; 
         p(1, nodes) .* p(2, nodes); 
         p(2, nodes).^2];
    
    areaOfElement = abs(det(P(1 : 3, 1 : 3))) / 2;
    
    % Three integration points within each triangle
    ip = eta_xi * p(:, nodes(1 : 3))';
 
    % 6 by 3 matrix with rows: [ones; x; y; x^2; xy; y^2]' (of three integration points)
    IPS = [ones(1, 3); ip'; ip(:, 1)'.^2; ip(:, 1)' .* ip(:, 2)'; ip(:, 2)'.^2];
    IPSPrimeX = [zeros(1, 3); ones(1, 3); zeros(1, 3); 2 * ip(:, 1)'; ip(:, 2)'; zeros(1, 3)];
    ISPPrimeY = [zeros(1, 3); zeros(1, 3); ones(1, 3); zeros(1, 3); ip(:, 1)'; 2 * ip(:, 2)'];
    % Concatenate
    IPSPrime =  [IPSPrimeX(:, 1) ISPPrimeY(:, 1) IPSPrimeX(:, 2) ISPPrimeY(:, 2) IPSPrimeX(:, 3) ISPPrimeY(:, 3)];    
    
    % Viscous terms
    PhiIPS = P \ IPS; 
    PhiIpsPrime = P \ IPSPrime;
    
    Ke = wOmega(1) * PhiIpsPrime(:, 1:2) *  PhiIpsPrime(:, 1:2)' * areaOfElement + ...
         wOmega(2) * PhiIpsPrime(:, 3:4) *  PhiIpsPrime(:, 3:4)' * areaOfElement + ...
         wOmega(3) * PhiIpsPrime(:, 5:6) * PhiIpsPrime(:, 5:6)' * areaOfElement;
     
    Me = wOmega(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * PhiIPS(:, 2) * PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * PhiIPS(:, 3) * PhiIPS(:, 3)' * areaOfElement;   
     
    %Convective Terms
    PhiDxIPS = P \ IPSPrimeX;
    PhiDyIPS = P \ ISPPrimeY;
    PhiDzIPS = P \ IPS; %no derivative because of the fourier transform
    
    Se1 = wOmega(1) * ubar(ip(1,:)) * PhiDxIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * ubar(ip(2,:)) * PhiDxIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * ubar(ip(3,:)) * PhiDxIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;

    Se2 = wOmega(1) * ubar(ip(1,:)) * PhiDyIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * ubar(ip(2,:)) * PhiDyIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * ubar(ip(3,:)) * PhiDyIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
     
    Te = wOmega(1) * ubar(ip(1,:)) * PhiDzIPS(:, 1) *  PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * ubar(ip(2,:)) * PhiDzIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * ubar(ip(3,:)) * PhiDzIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
       
    
    %Pressure component
    %phis are quadratic, psis are linear
    psiPres= P(1:3,1:3)\IPS(1:3,1:3);
    %PhiDxIPS = P \ IPSPrimeX;
    %PhiDyIPS = P \ ISPPrimeY;
    %PhiDzIPS = P \ IPS; %no derivative because of the fourier transform
    
    Be1 = wOmega(1) * PhiDxIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
         wOmega(2) * PhiDxIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
         wOmega(3) * PhiDxIPS(:, 3) * psiPres(:, 3)' * areaOfElement;

    Be2 = wOmega(1) * PhiDyIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
         wOmega(2) * PhiDyIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
         wOmega(3) * PhiDyIPS(:, 3) * psiPres(:, 3)' * areaOfElement;
     
    Be3 = wOmega(1) * PhiDzIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
         wOmega(2) * PhiDzIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
         wOmega(3) * PhiDzIPS(:, 3) * psiPres(:, 3)' * areaOfElement;
     
     
      nodes12 = [nodes, nodes+numberOfNodes.new];
      nodesOld = nodes(1:3);
      K(nodes, nodes) = K(nodes, nodes) + Ke;
      M(nodes, nodes) = M(nodes, nodes) + Me;
      
      S1(nodes, nodes) = S1(nodes, nodes) + Se1;
      S2(nodes, nodes) = S2(nodes, nodes) + Se2;
      T(nodes, nodes) = T(nodes, nodes) + Te;
      
      B(nodes12, nodesOld) = B(nodes12, nodesOld) + [Be1; Be2];
      B3(nodes, nodesOld) = B3(nodes, nodesOld)+Be3;
      
      
      
    nodes123=[nodes, nodes+numberOfNodes.new, nodes+2*numberOfNodes.new];
    for i=1:maxWaveNum
        waveNum=waveNumbers(i);
        Fe = wOmega(1) * PhiIPS(:, 1) * f(ip(1, :), waveNum) * areaOfElement + ...
         wOmega(2) * PhiIPS(:, 2) * f(ip(2, :), waveNum) * areaOfElement + ...
         wOmega(3) * PhiIPS(:, 3) * f(ip(3, :), waveNum) * areaOfElement; 
        nodesi=nodes123+(i-1)*(wnl);
        F(nodesi) = F(nodesi) + reshape(Fe, 3*length(nodes), 1);   
    end
end

% Putting together the block function
Bs=[sparse(B);sparse(numberOfNodes.new, numberOfNodes.old)];
B3s=[sparse(2*numberOfNodes.new, numberOfNodes.old); sparse(B3)];
M=sparse(M);
K=sparse(K);
Btot = @(waveNum) Bs + 1i*(2*pi)^0.5*waveNum*B3s;
Visc =@(waveNum) -K+2*pi*waveNum.^2*M;
momentumEQNs = @(waveNum) [Visc(waveNum), sparse(numberOfNodes.new,numberOfNodes.new), 1i*(2*pi)^0.5*waveNum*T;
                           sparse(numberOfNodes.new,numberOfNodes.new), Visc(waveNum), 1i*(2*pi)^0.5*waveNum*T;
                           S1, S2, Visc(waveNum)+ 1i*(2*pi)^0.5*waveNum*T];
    
Awn = @(waveNum) [momentumEQNs(waveNum), Btot(waveNum); Btot(waveNum)', sparse(numberOfNodes.old,numberOfNodes.old)];

Awns=cell(1,length(waveNumbers));
PCG = cell(1,length(waveNumbers));
for i=1:length(waveNumbers)
    Awns(i) = {Awn(waveNumbers(i))};
end
%rank one update to increase rank
u = [zeros(3*numberOfNodes.new,1);ones(numberOfNodes.old,1)];
u = sparse(u);
Awns{maxWaveNum/2} = Awns{maxWaveNum/2}+1*mean(mean(abs(Bs)))*(u*u');


% Dirichlet boundary
Dirichlet123 = [Dirichlet;Dirichlet+numberOfNodes.new;Dirichlet+2*numberOfNodes.new];
for i=1:maxWaveNum
    Awns{i}(Dirichlet123, :) = 0;
    Awns{i}(Dirichlet123, Dirichlet123) = eye(numel(Dirichlet123));
    
    F(Dirichlet123+(i-1)*wnl) = Fb((1:numel(Dirichlet123))+(i-1)*numel(Dirichlet123));
end

%A = blkdiag(Awns{:});
clear K M Bs B3s
U=cell(maxWaveNum);
for i=1:maxWaveNum
    U{i} = Awns{i}\F((1:wnl)+(i-1)*wnl);
end
    
    


tol=10^(-10);
maxit=100;
%U=gmres(A,F,30,tol,maxit, PCG);
Uwn=cell(maxWaveNum, 4); % (waveNumbers, u1 u2 u3  p) 
for i=1:maxWaveNum
   Uwn{i,1} = U{i}((1:numberOfNodes.old)); 
   Uwn{i,2} = U{i}((1:numberOfNodes.old)+numberOfNodes.new);
   Uwn{i,3} = U{i}((1:numberOfNodes.old)+2*numberOfNodes.new);
   Uwn{i,4} = U{i}((1:numberOfNodes.old)+3*numberOfNodes.new);
end
%for i=1:maxWaveNum
%   Uwn{i,1} = U((1:numberOfNodes.old)+(i-1)*wnl); 
%   Uwn{i,2} = U((1:numberOfNodes.old)+(i-1)*wnl+numberOfNodes.new);
%   Uwn{i,3} = U((1:numberOfNodes.old)+(i-1)*wnl+2*numberOfNodes.new);
%   Uwn{i,4} = U((1:numberOfNodes.old)+(i-1)*wnl+3*numberOfNodes.new);
%end
