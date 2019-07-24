function [velocity,pOld,tOld] = velocitySolve(p,t,epsilon,waveNumbers,xp)
%% Calculate stuff used in all solves.
maxWaveNum=length(waveNumbers);
numPart = size(xp, 2);
b=boundedges(p,t);
Dirichlet_e=b;
p=p';
fb = @(x, k) [0*x(:,2), 0*x(:,1), 0*x(:,1)+0*k];

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
RHSfunc = cell(numPart, 1);
us = zeros(numPart,1);
gammax = zeros(1, numPart);
gammay = zeros(1, numPart,1);
bgFlow = 1-p(1, :).^2-p(2, :).^2;
particleNodes = cell(numPart,1);
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

N = @(x,k,xpf) 100*epsilon^2/((epsilon^2*pi*2)^1.5)*exp(-1/2*(((x(:,1)-xpf(1,:)).^2+(x(:,2)-xpf(2,:)).^2)/epsilon^2+epsilon^2*k.^2));
f=@(x,xpf,k,gx,gy) [-gx.*1i*k.*N(x, k, xpf),...
    -gy.*1i*k.*N(x,k, xpf), 1/epsilon^2*(gx.*(x(:,1)-xpf(1,:)).*N(x,k,xpf)+gy.*(x(:,2)-xpf(2,:)).*N(x,k,xpf))];



%
% N = @(x,k) epsilon/((epsilon^2*pi*2)^1.5)*exp(-1/2*(((x(:,1)-xp(1)).^2+(x(:,2)-xp(2)).^2)/epsilon^2+epsilon^2*k^2));
% gammax = bgFlowIPS * PhiDxIPS;
% gammay = bgFlowIPS * PhiDyIPS;
% f=@(x,k) [-gammax*N(x, k),...
%           -gammay*N(x,k), -(gammax*N(x,k)+gammay*N(x,k))];



% Initialisation of K, M, B, S1, S2, T, F
K = sparse(numberOfNodes.new,numberOfNodes.new);
M = sparse(numberOfNodes.new,numberOfNodes.new);
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
Btot = @(waveNum) Bs + 1i*waveNum*B3s;
Visc =@(waveNum) K+waveNum.^2*M;



% Create Preconditioner
a1 = 3*numberOfNodes.new;
M(Dirichlet, :) = 0;
K(Dirichlet, :) = 0;
K(Dirichlet, Dirichlet) = -eye(numel(Dirichlet));
Visc = @(waveNum) K + waveNum.^2*M;
Fhat = @(waveNum) blkdiag(-Visc(waveNum), -Visc(waveNum), -Visc(waveNum));

%% Do individual solves for particles
momentumEQNs = @(waveNum, us) [K + (waveNum.^2-1i*waveNum*us)*M + 1i*waveNum*T, sparse(numberOfNodes.new, 2*numberOfNodes.new);
    sparse(numberOfNodes.new,numberOfNodes.new), K + (waveNum.^2-1i*waveNum*us)*M + 1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
    S1, S2, K + (waveNum.^2-1i*waveNum*us)*M + 1i*waveNum*T];

Awn = @(waveNum, us) [-momentumEQNs(waveNum, us), Btot(waveNum); Btot(waveNum)', sparse(numberOfNodes.old,numberOfNodes.old)];

velocity = zeros(numPart, 2);
%AwnPool = parallel.pool.Constant(@() Awn);

UhatNodes1 = zeros(6, maxWaveNum, numPart);
UhatNodes2 = zeros(6, maxWaveNum, numPart);
UhatNodes3 = zeros(6, maxWaveNum, numPart);
parfor waveIndex = 1:maxWaveNum
    
    
    Fb=zeros(length(Dirichlet)*numPart,1);
    for i=1:numPart
        %this part is dumb but doesn't effect anything because it's all
        %zeros. I keep it because it might be useful if I want to implement
        %non zero boundary conditions.
        waveNum=waveNumbers(i);
        Fb((1:3*length(Dirichlet))+(i-1)*3*length(Dirichlet))=reshape(fb(p(:,Dirichlet)',waveNum),3*length(Dirichlet),1);
    end
    
    F = zeros(wnl, numPart);
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
        
        
        % right hand side
        nodes123 = [nodes, nodes+numberOfNodes.new, nodes+2*numberOfNodes.new];
        
        %linearize with respect to multiple xp
        PhiIPS = P \ IPS;
        Fe = wOmega(1) * PhiIPS(:, 1) * f(ip(1, :), xp, waveNumbers(waveIndex), gammax, gammay ) * areaOfElement + ...
            wOmega(2) * PhiIPS(:, 2) * f(ip(2, :), xp, waveNumbers(waveIndex), gammax, gammay) * areaOfElement + ...
            wOmega(3) * PhiIPS(:, 3) * f(ip(3, :), xp, waveNumbers(waveIndex), gammax, gammay) * areaOfElement;
        
        F(nodes123, :) = F(nodes123, :) + [Fe(:,1:end/3); Fe(:,(1:end/3) + end/3); Fe(:,(1:end/3) + 2*end/3)];
    end
    
    %     Awns=cell(1,length(waveNumbers));
    %     for i=1:length(waveNumbers)
    %         Awns(i) = {Awn(waveNumbers(i))};
    %     end
    %rank one update to increase rank
    u = [zeros(3*numberOfNodes.new,1); ones(numberOfNodes.old,1)];
    u = sparse(u);
    %     Awns{1} = Awns{1}+100*mean(mean(abs(Bs)))*(u*u');
    
    
    % Dirichlet boundary
    Dirichlet123 = [Dirichlet; Dirichlet+numberOfNodes.new;Dirichlet+2*numberOfNodes.new];
    %     for i=1:maxWaveNum
    %         %         Awns{i}(Dirichlet123, :) = 0;
    %         %         Awns{i}(Dirichlet123, Dirichlet123) = eye(numel(Dirichlet123));
    %         %
    %         F(Dirichlet123, i) = Fb((1:numel(Dirichlet123))+(i-1)*numel(Dirichlet123));
    %     end
    
    % for i=1:maxWaveNum
    %     U{i} = Awns{i}\F((1:wnl)+(i-1)*wnl);
    % end
    [FhatL, FhatU,FhatP,FhatQ] = lu(Fhat(waveNumbers(waveIndex)));
    schurCompi = 100*mean(mean(abs(Bs)))*eq(waveIndex,1)*ones(numberOfNodes.old)+Mp;
    [schurL, schurU] = lu(schurCompi);
    tol = 10^(-6);
    maxit = 30;
    buttshit = zeros(6, 3, numPart);
    for i = 1:numPart
        A = Awn(waveNumbers(waveIndex), us(i));
        A(Dirichlet123, :) = 0;
        A(Dirichlet123, Dirichlet123) = eye(numel(Dirichlet123));
        if i ==1
            A = A+100*mean(mean(abs(Bs)))*(u*u');
        end
        
        Uhat = gmres(A, F(:, i),15, tol,maxit, @(x) PCbackSolve(x,A(1:a1, (a1+1):end),FhatL,FhatU, FhatP, FhatQ, schurL, schurU, a1));
        buttshit(:,1, i) = Uhat(particleNodes{i});
        buttshit(:,2, i) = Uhat(particleNodes{i} + numberOfNodes.new);
        buttshit(:,3, i) = Uhat(particleNodes{i} + 2*numberOfNodes.new);
    end
    UhatNodes1(:, waveIndex, :) = buttshit(:, 1, :);
    UhatNodes2(:, waveIndex, :) = buttshit(:, 2, :);
    
    %     Uwn=cell(maxWaveNum, 4); % (waveNumbers, u1 u2 u3  p)
    %     for i=1:maxWaveNum
    %         Uwn{i,1} = U{i}((1:numberOfNodes.new));
    %         Uwn{i,2} = U{i}((1:numberOfNodes.new)+numberOfNodes.new);
    %         Uwn{i,3} = U{i}((1:numberOfNodes.new)+2*numberOfNodes.new);
    %         Uwn{i,4} = U{i}((1:numberOfNodes.old)+3*numberOfNodes.new);
    %     end

end

%inverse fourier transform and evaluate velocity

Unodes1 = ifft(UhatNodes1, maxWaveNum, 2);
Unodes2 = ifft(UhatNodes2, maxWaveNum, 2);
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
    
    velocity(i,1) = Unodes1(:,1,i)' * PhiIPS;
    velocity(i,2) = Unodes2(:,1,i)' * PhiIPS;
end
