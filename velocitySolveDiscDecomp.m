function [velocity,pOld,tOld] = velocitySolveDiscDecomp(p,t,epsilon,waveNumbers,xp,bgFlow)
%% Calculate stuff used in all solves.
maxWaveNum=length(waveNumbers);
numPart = size(xp, 2);
b=boundedges(p,t);
Dirichlet_e=b;
p=p';

pOld = p;
tOld = t;

% Quadratic elements contain six nodes per triangle: add three nodes at the middle of the edges of each triangle
numberOfNodes.old = size(p, 2);
numberOfElements = size(t, 1);
S = sparse(numberOfNodes.old, numberOfNodes.old);
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

% creat f for each particle
% Frame things in terms of xp
us = zeros(numPart,1);
gammax = zeros(1, numPart);
gammay = zeros(1, numPart,1);

particleNodes = cell(numPart,1);
uMax = max(bgFlow);
bgFlow = bgFlow/uMax;
for part =1:numPart
    xpi = xp(:, part);
    xpn = whatTriangleIsThisPointIn(p,t,xpi);
    nodes = t(xpn, :);
    particleNodes{part} = nodes;
    % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
    P = [ones(1, 6);
        p(:, nodes);
        p(1, nodes).^2;
        p(1, nodes) .* p(2, nodes);
        p(2, nodes).^2];
    IPS = [1; xpi; xpi(1)^2; xpi(1) * xpi(2); xpi(2)^2];
    IPSPrimeX = [zeros(1, 1); ones(1, 1); zeros(1, 1); 2 * xpi(1); xpi(2); zeros(1, 1)];
    ISPPrimeY = [zeros(1, 1); zeros(1, 1); ones(1, 1); zeros(1, 1); xpi(1); 2 * xpi(2)];
    
    PhiIPS = P \ IPS;
    PhiDxIPS = P \ IPSPrimeX;
    PhiDyIPS = P \ ISPPrimeY;
    
    bgFlowIPS = bgFlow(nodes);
    us(part) = bgFlowIPS * PhiIPS;
    gammax(part) = bgFlowIPS * PhiDxIPS;
    gammay(part) = bgFlowIPS * PhiDyIPS;
end
ubar = bgFlow;

% Initialisation of K, M, B, S1, S2, T, F0
K = sparse(numberOfNodes.new,numberOfNodes.new);
M = sparse(numberOfNodes.new,numberOfNodes.new);
Mx = sparse(numberOfNodes.new,numberOfNodes.new);
My = sparse(numberOfNodes.new,numberOfNodes.new);
Mp = sparse(numberOfNodes.old,numberOfNodes.old);
B = sparse(2*numberOfNodes.new,numberOfNodes.old);
B3 = sparse(numberOfNodes.new,numberOfNodes.old);
S1 = sparse(numberOfNodes.new, numberOfNodes.new);
S2 = sparse(numberOfNodes.new, numberOfNodes.new);
T = sparse(numberOfNodes.new, numberOfNodes.new);

% Gaussian quadrature points & weights
eta_xi = [2/3 1/6 1/6;
    1/6 2/3 1/6;
    1/6 1/6 2/3];
wOmega = [1/3 1/3 1/3];
%% The for loop
% Assembly of K, M, and B
% Assembly of generic S1, S2, and T
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
    
    Mxe = p(1, nodes(1)) * wOmega(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        p(1, nodes(2)) * wOmega(2) * PhiIPS(:, 2) * PhiIPS(:, 2)' * areaOfElement + ...
        p(1, nodes(3)) * wOmega(3) * PhiIPS(:, 3) * PhiIPS(:, 3)' * areaOfElement;
    
    Mye = p(2, nodes(1)) * wOmega(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        p(2, nodes(2)) * wOmega(2) * PhiIPS(:, 2) * PhiIPS(:, 2)' * areaOfElement + ...
        p(2, nodes(3)) * wOmega(3) * PhiIPS(:, 3) * PhiIPS(:, 3)' * areaOfElement;
    
    %Pressure component
    psiPres = P(1:3,1:3)\IPS(1:3,1:3);
    PhiDxIPS = P \ IPSPrimeX;
    PhiDyIPS = P \ ISPPrimeY;
    PhiDzIPS = P \ IPS; %no derivative because of the fourier transform
    
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
    Mx(nodes, nodes) = Mx(nodes, nodes) + Mxe;
    My(nodes, nodes) = My(nodes, nodes) + Mye;
    
    Mp(nodesOld, nodesOld) = Mp(nodesOld, nodesOld) + Mpe;
    
    B(nodes12, nodesOld) = B(nodes12, nodesOld) + [Be1; Be2];
    B3(nodes, nodesOld) = B3(nodes, nodesOld)+Be3;
    
    
    % Calculate S1, S2, and T
    
    
    
    
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
    
    % Viscous terms
    PhiIPS = P \ IPS;
    
    %Convective Terms
    PhiDxIPS = P \ IPSPrimeX;
    PhiDyIPS = P \ ISPPrimeY;
    PhiDzIPS = P \ IPS; %no derivative because of the fourier transform
    
    
    ubarx_ip = ubar(nodes)*PhiDxIPS;
    ubary_ip = ubar(nodes)*PhiDyIPS;
    ubar_ip = ubar(nodes)*PhiIPS;
    
    
    Se1 = wOmega(1) * ubarx_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubarx_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubarx_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
    Se2 = wOmega(1) * ubary_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubary_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubary_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
    Te = wOmega(1) * ubar_ip(1) * PhiDzIPS(:, 1) *  PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubar_ip(2) * PhiDzIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubar_ip(3) * PhiDzIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
    S1(nodes, nodes) = S1(nodes, nodes) + Se1;
    S2(nodes, nodes) = S2(nodes, nodes) + Se2;
    T(nodes, nodes) = T(nodes, nodes) + Te;
    
end

% Putting together the block function
Bs=[sparse(B);sparse(numberOfNodes.new, numberOfNodes.old)];
B3s=[sparse(2*numberOfNodes.new, numberOfNodes.old); sparse(B3)];
M=sparse(M);
K=sparse(K);
a1 = 3*numberOfNodes.new;
M(Dirichlet, :) = 0;
T(Dirichlet, :) = 0;
S1(Dirichlet, :) = 0;
S2(Dirichlet, :) = 0;
K(Dirichlet, :) = 0;
K(Dirichlet, Dirichlet) = -eye(numel(Dirichlet));




%% Do individual solves for particles

velocity = zeros(numPart, 2);

UhatNodes1 = zeros(6, maxWaveNum, numPart);
UhatNodes2 = zeros(6, maxWaveNum, numPart);
UhatNodes3 = zeros(6, maxWaveNum, numPart);

Btot = @(waveNum) Bs + 1i*waveNum*B3s;
momentumEQNs = @(waveNum) [K + (waveNum.^2)*M - 1i*waveNum*T, sparse(numberOfNodes.new, 2*numberOfNodes.new);
    sparse(numberOfNodes.new,numberOfNodes.new), K + (waveNum.^2)*M - 1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
    S1, S2, K + (waveNum.^2)*M - 1i*waveNum*T];

% Create Preconditioner
Visc = @(waveNum) K + waveNum.^2*M;
Fhat = @(waveNum) blkdiag(-Visc(waveNum), -Visc(waveNum), -Visc(waveNum));

AwnBase = @(waveNum) [-momentumEQNs(waveNum), -Btot(waveNum); -Btot(waveNum)', sparse(numberOfNodes.old,numberOfNodes.old)];
uStressF = cell(numPart);
discF = cell(numPart);

for i = 1:numPart
uStressF{i} = regStressletF(p',xp(:,i),gammax(i),gammay(i),epsilon, 2*maxWaveNum+2);
uStressF{i} = uStressF{i}(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/4*2);
discF{i} = uDiscF(p', xp(1,i), xp(2,i), gammax(i), gammay(i), 2*maxWaveNum+2);
discF{i} = 0*discF{i}(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/4*2);
end

RHSop = @(waveNum) [-1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new*2);
                    sparse(numberOfNodes.new,numberOfNodes.new), -1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
                    S1, S2 , -1i*waveNum*T];
RHSstress = @(waveNum,i) [+1i*waveNum*(us(i)-gammax(i)*xp(1))*M, sparse(numberOfNodes.new,numberOfNodes.new*2);
                    sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*(us(i)-gammay(i)*xp(2))*M, sparse(numberOfNodes.new,numberOfNodes.new);
                    sparse(numberOfNodes.new, 2*numberOfNodes.new) , 1i*waveNum*us(i)*M];

RHSstress2 = @(waveNum,i) [+1i*waveNum*gammax(i)*Mx, sparse(numberOfNodes.new,numberOfNodes.new*2);
                    sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*gammay(i)*My, sparse(numberOfNodes.new,numberOfNodes.new);
                    -gammax(i)*M, -gammay(i)*M, sparse(numberOfNodes.new,numberOfNodes.new)];              
cellforblkdiag = {M,M,M}; %apparently this is the best way to do this                
for waveIndex = 1:maxWaveNum
    
    
    %rank one update to increase rank
    u = [zeros(3*numberOfNodes.new,1); ones(numberOfNodes.old,1)];
    u = sparse(u);
   
    
    % Dirichlet boundary
    Dirichlet123 = [Dirichlet; Dirichlet+numberOfNodes.new;Dirichlet+2*numberOfNodes.new];
  
    
    [FhatL, FhatU,FhatP,FhatQ] = lu(Fhat(waveNumbers(waveIndex)));
    schurCompi = 100*mean(mean(abs(Bs)))*eq(waveIndex,1)*ones(numberOfNodes.old)+Mp;
    [schurL, schurU] = lu(schurCompi);
    tol = 10^(-4);
    maxit = 30;
    buttshit = zeros(6, 3, numPart); %you need to create this because of how matlab does slicing
    RHSbase = RHSop(waveNumbers(waveIndex));
    F = zeros(wnl,1);
    ABase = AwnBase(waveNumbers(waveIndex));
    A=ABase;
    A(Dirichlet123, :) = 0;
    for i = 1:numPart
        A(1:a1/3,1:a1/3) = ABase(1:a1/3,1:a1/3)-1i*waveNumbers(waveIndex)*us(i)*M;
        A((a1/3+1):(2*a1/3),(a1/3+1):(2*a1/3)) = A(1:a1/3,1:a1/3);
        A((2*a1/3+1):(a1),(2*a1/3+1):(a1)) = A(1:a1/3,1:a1/3);
         F(:,i) = [RHSbase*(uStressF{i}(:,waveIndex)-discF{i}(:,waveIndex))+RHSstress(waveNumbers(waveIndex),i)*(uStressF{i}(:,waveIndex))...
                  + RHSstress2(waveNumbers(waveIndex),i)*uStressF{i}(:,waveIndex)...
                  + 1i*waveNumbers(waveIndex)*us(i)*blkdiag(cellforblkdiag{:})*discF{i}(:,waveIndex);
                   sparse(numberOfNodes.old,1)]; 
         F(Dirichlet123, i) = -uStressF{i}(Dirichlet123,waveIndex) + discF{i}(Dirichlet123,waveIndex);
        %A(Dirichlet123, :) = 0;
        %A(Dirichlet123, Dirichlet123) = eye(numel(Dirichlet123));
        if waveIndex == 1
            A = A + 100*mean(mean(abs(Bs)))*(u*u');
        end
        
        Uhat = gmres(A, F(:, i),15, tol,maxit, @(x) PCbackSolve(x,A(1:a1, (a1+1):end),FhatL,FhatU, FhatP, FhatQ, schurL, schurU, a1));
    
        buttshit(:,1, i) = Uhat(particleNodes{i});
        buttshit(:,2, i) = Uhat(particleNodes{i} + numberOfNodes.new);
        buttshit(:,3, i) = Uhat(particleNodes{i} + 2*numberOfNodes.new);
    end
    UhatNodes1(:, waveIndex, :) = buttshit(:, 1, :);
    UhatNodes2(:, waveIndex, :) = buttshit(:, 2, :);
    UhatNodes3(:, waveIndex, :) = buttshit(:, 3, :);

end

%inverse fourier transform and evaluate velocity

Unodes1 = sum(real(UhatNodes1), 2)/maxWaveNum;
Unodes2 = sum(real(UhatNodes2), 2)/maxWaveNum;
for i = 1:numPart
    nodes = particleNodes{i};
    xpi = xp(:,i);
    P = [ones(1, 6);
        p(:, nodes);
        p(1, nodes).^2;
        p(1, nodes) .* p(2, nodes);
        p(2, nodes).^2];
    IPS = [1; xpi; xpi(1)^2; xpi(1) * xpi(2); xpi(2)^2];
    
    PhiIPS = P \ IPS;
    
   velocity(i,1) = Unodes1(:,i)' * PhiIPS;
   velocity(i,2) = Unodes2(:,i)' * PhiIPS;
end