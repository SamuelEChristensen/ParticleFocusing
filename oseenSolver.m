function [Uwn,pOld,tOld] = oseenSolver(p,t,f,fb,waveNumbers,xp)

maxWaveNum=length(waveNumbers);
b=boundedges(p,t);
Dirichlet_e=b;
p=p';

 ubar =@(x) xp(1)^2+xp(2)^2 - x(:,1).^2-x(:,2).^2 ;

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


% figure out column swaps to make matrix more diagonal
% R = zeros(numberOfNodes.new - numberOfNodes.old,1);
% count = 0;
% for i = (numberOfNodes.old + 1):numberOfNodes.new
%     count = count + 1;
%     [a, ~] = find(S == i);
%     R(count) = mean(a);
% end
% R = [(1:numberOfNodes.old)'; R];
% [~, I] = sort(R);
I = 1:numberOfNodes.new;

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
Mp = sparse(numberOfNodes.old,numberOfNodes.old);
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
    
    
    ubarx_ip = ubar(p(:,nodes)')'*PhiDxIPS;
    ubary_ip = ubar(p(:,nodes)')'*PhiDyIPS;

    
   Se1 = wOmega(1) * ubarx_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * ubarx_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * ubarx_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;

    Se2 = wOmega(1) * ubary_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
         wOmega(2) * ubary_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
         wOmega(3) * ubary_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;     

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
     
    Mpe = wOmega(1) * psiPres(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
         wOmega(2) * psiPres(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
         wOmega(3) * psiPres(:, 3) * psiPres(:, 3)' * areaOfElement;
     
     
      nodes12 = [nodes, nodes+numberOfNodes.new];
      nodesOld = nodes(1:3);
      K(nodes, nodes) = K(nodes, nodes) + Ke;
      M(nodes, nodes) = M(nodes, nodes) + Me;
      
      Mp(nodesOld, nodesOld) = Mp(nodesOld, nodesOld) + Mpe;
      
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
Btot = @(waveNum) Bs + 1i*waveNum*B3s;
Visc =@(waveNum) K+waveNum.^2*M;
momentumEQNs = @(waveNum) [Visc(waveNum) + 1i*waveNum*T, sparse(numberOfNodes.new, 2*numberOfNodes.new);
                           sparse(numberOfNodes.new,numberOfNodes.new), Visc(waveNum) + 1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
                           S1, S2, Visc(waveNum) + 1i*waveNum*T];
    
Awn = @(waveNum) [-momentumEQNs(waveNum), Btot(waveNum); Btot(waveNum)', sparse(numberOfNodes.old,numberOfNodes.old)];

Awns=cell(1,length(waveNumbers));
for i=1:length(waveNumbers)
    Awns(i) = {Awn(waveNumbers(i))};
end


%rank one update to increase rank
u = [zeros(3*numberOfNodes.new,1); ones(numberOfNodes.old,1)];
u = sparse(u);
Awns{1} = Awns{1}+100*mean(mean(abs(Bs)))*(u*u');

% Dirichlet boundary and diagonlizing
%Dirichlet123s = [I(Dirichlet); I(Dirichlet)+numberOfNodes.new;I(Dirichlet)+2*numberOfNodes.new];
Dirichlet123 = [Dirichlet; Dirichlet+numberOfNodes.new;Dirichlet+2*numberOfNodes.new];
switchbackm = speye(length(I));
switchbackm = switchbackm(I,:);
switchbackp = speye(numberOfNodes.old);
switchback = blkdiag(switchbackm,switchbackm,switchbackm,switchbackp);

for i=1:maxWaveNum
     Awns{i}(Dirichlet123, :) = 0;
     Awns{i}(Dirichlet123, Dirichlet123) = eye(numel(Dirichlet123));
     Awns{i} = switchback*Awns{i}*switchback';
    F(Dirichlet123+(i-1)*wnl) = Fb((1:numel(Dirichlet123))+(i-1)*numel(Dirichlet123));
end

% Create Preconditioner
a1 = 3*numberOfNodes.new;
M(Dirichlet, :) = 0;
K(Dirichlet, :) = 0;
K(Dirichlet, Dirichlet) = -eye(numel(Dirichlet));
T(Dirichlet, :) = 0;
K = switchbackm*K*switchbackm';
M = switchbackm*M*switchbackm';
T = switchbackm*T*switchbackm';
Visc = @(waveNum) K+waveNum.^2*M + 0*1i*waveNum*T;
Fhat = @(waveNum) blkdiag(-Visc(waveNum),-Visc(waveNum),-Visc(waveNum));



%schurComp = @(index) Awns{index}((a1+1):end,1:a1)*(Awns{index}(1:a1,1:a1)\Awns{index}(1:a1, (1+a1):end));

%PCR = @(index) [Awns{index}(1:a1, 1:a1), Awns{index}(1:a1, (a1+1):end); zeros(numberOfNodes.old, a1), 100*mean(mean(abs(Bs)))*eq(index,1)*ones(numberOfNodes.old)+schurComp(index)];
%PCL = @(index) [Awns{index}(1:a1,1:a1), zeros(a1, numberOfNodes.old);Awns{index}((a1+1):end, 1:a1), 100*mean(mean(abs(Bs)))*eq(index,1)*ones(numberOfNodes.old)+schurComp(index)];
schurCompsL = cell(maxWaveNum);
schurCompsU = cell(maxWaveNum);
FhatL = cell(maxWaveNum);
FhatU = cell(maxWaveNum);
FhatP = cell(maxWaveNum);
FhatQ = cell(maxWaveNum);
schurComp = @(index,L,U,P,Q) Awns{index}((a1+1):end,1:a1)*(Q*(U\(L\(P*Awns{index}(1:a1, (1+a1):end)))));

for i = 1:maxWaveNum
   [L1,U1,P1,Q1] = lu(Fhat(waveNumbers(i)));
   FhatL(i) = {L1};
   FhatU(i) = {U1};
   FhatP(i) = {P1};
   FhatQ(i) = {Q1};
   schurCompi = 100*mean(mean(abs(Bs)))*eq(i,1)*ones(numberOfNodes.old)+Mp;%100*mean(mean(abs(Bs)))*eq(i,1)*ones(numberOfNodes.old)-schurComp(i,L1,U1,P1,Q1);
   [L2,U2] = lu(schurCompi);
   schurCompsL{i} = L2;
   schurCompsU{i} = U2;
end

clear K M Bs B3s
U=cell(maxWaveNum);

% for i=1:maxWaveNum
% %      Awns{i} = switchback*Awns{i}*switchback';
%       U{i} = switchback'*(Awns{i}\(switchback*F((1:wnl)+(i-1)*wnl)));
% %    U{i} = Awns{i}\F((1:wnl)+(i-1)*wnl);
% end
    
    


tol = 10^(-9);
maxit = 30;
tic
for i = 1:maxWaveNum
U{i} = gmres(Awns{i},switchback*F((1:wnl)+(i-1)*wnl),15, tol,maxit, @(x) PCbackSolve(x,Awns{i}(1:a1, (a1+1):end),FhatL{i}, FhatU{i}, FhatP{i}, FhatQ{i}, schurCompsL{i}, schurCompsU{i},a1));
end

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
end
