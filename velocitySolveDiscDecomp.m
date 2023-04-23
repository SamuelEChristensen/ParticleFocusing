function [velocity,p,t] = velocitySolveDiscDecomp(p,t,epsilon,waveNumbers,xp,bgFlow,L,nParam, D, Re)
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
S = zeros(numberOfNodes.old, numberOfNodes.old);%sparse(numberOfNodes.old, numberOfNodes.old);
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
clear S
% creat f for each particle
% Frame things in terms of xp
us = zeros(numPart,1);
gammax = zeros(1, numPart);
gammay = zeros(1, numPart,1);

particleNodes = cell(numPart,1);
nearbyTris = cell(numPart,1);
uMax = max(bgFlow);
bgFlow = Re * bgFlow/uMax;
for part =1:numPart
    xpi = xp(:, part);
    xpn = whatTriangleIsThisPointIn(p,t,xpi);
    nearbyTris{part}=[xpn];
    for i =1:D
        nearbyTris{part} = nearbyTriangles(t,t(nearbyTris{part},:));  %can be iterated for more surrounding layers
        if i== (D-1)
            intNearbyTris{part} = nearbyTris{part};
        end
        if D ==1
            intNearbyTris{part} = xpn;
        end
        
    end
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
    gammax(part) = bgFlowIPS * PhiDxIPS/Re;
    gammay(part) = bgFlowIPS * PhiDyIPS/Re;
end
ubar = bgFlow;
%% Dealing with subtriangles
subtp = cell(numPart,1);
subtt = cell(numPart,1);
Mcell = cell(numPart,1);
Mxcell = cell(numPart,1);
Mycell = cell(numPart,1);
S1cell = cell(numPart,1);
S2cell = cell(numPart,1);
Tcell = cell(numPart,1);
% Gaussian quadrature points & weights
eta_xi = [2/3 1/6 1/6;
    1/6 2/3 1/6;
    1/6 1/6 2/3];
wOmega = [1/3 1/3 1/3];

% Initialisation of K, M, S1, S2, T, F0 for sub triangles

tic
for part = 1:numPart
    tCell = cell(numel(nearbyTris{part}), 1);
    pCell = cell(numel(nearbyTris{part}), 1);
    Msubcell = cell(numel(nearbyTris{part}), 1);
    Mxsubcell = cell(numel(nearbyTris{part}), 1);
    Mysubcell = cell(numel(nearbyTris{part}), 1);
    Tsubcell = cell(numel(nearbyTris{part}), 1);
    S1subcell = cell(numel(nearbyTris{part}), 1);
    S2subcell = cell(numel(nearbyTris{part}), 1);
    
    for e = 1:numel(nearbyTris{part})
        
        [subp,subt] = subTriangleMaker(p(:,t(nearbyTris{part}(e),1:3))');
        allp=subp;
        allt=subt;
        for j = 1:nParam
            allpi=[];
            allti=[];
            for i = 1:length(allt)  %iterate this part more for more subtriangles
                [ssp,sst] = subTriangleMaker(allp(allt(i,:),:));
                allpi =[allpi; ssp];
                allti=[allti; sst+6*(i-1)];
            end
            allp = [allpi];
            allt = [allti];
        end
        supNodes = t(nearbyTris{part}(e), :);
        
        % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
        supP = [ones(1, 6);
            p(:, supNodes);
            p(1, supNodes).^2;
            p(1, supNodes) .* p(2, supNodes);
            p(2, supNodes).^2];
        
        allp=allp';
        F = size(allt,1);
        E = size(allp,2);
        [allp,allt] = makeMidPoints(allp,allt,F,E);
        F = size(allt,1);
        E = size(allp,2);
        Ms = zeros(E,6);
        Mxs = zeros(E,6);
        Mys = zeros(E,6);
        S1s = zeros(E,6);
        S2s = zeros(E,6);
        Ts = zeros(E,6);
        for sube = 1:F
            infNodes = allt(sube,:);
            infP = [ones(1, 6);
                allp(:, infNodes);
                allp(1, infNodes).^2;
                allp(1, infNodes) .* allp(2, infNodes);
                allp(2, infNodes).^2];
            areaOfElement = abs(det(infP(1 : 3, 1 : 3))) / 2;
            % Three integration points within each triangle
            ip = eta_xi * allp(:,infNodes(1:3))';
            
            % 6 by 3 matrix with rows: [ones; x; y; x^2; xy; y^2]' (of three integration points)
            IPS = [ones(1, 3); ip'; ip(:, 1)'.^2; ip(:, 1)' .* ip(:, 2)'; ip(:, 2)'.^2];
            IPSPrimeX = [zeros(1, 3); ones(1, 3); zeros(1, 3); 2 * ip(:, 1)'; ip(:, 2)'; zeros(1, 3)];
            ISPPrimeY = [zeros(1, 3); zeros(1, 3); ones(1, 3); zeros(1, 3); ip(:, 1)'; 2 * ip(:, 2)'];
            %            % Concatenate
            IPSPrime =  [IPSPrimeX(:, 1) ISPPrimeY(:, 1) IPSPrimeX(:, 2) ISPPrimeY(:, 2) IPSPrimeX(:, 3) ISPPrimeY(:, 3)];
            
            
            supPhiIPS = supP\IPS;
            %supPhiIpsPrime = supP \ IPSPrime;
            supPhiDxIPS = supP \ IPSPrimeX;
            supPhiDyIPS = supP \ ISPPrimeY;
            
            infPhiIPS = infP\IPS;
            %infPhiIpsPrime = infP \ IPSPrime;
            %infPhiDxIPS = infP \ IPSPrimeX;
            %infPhiDyIPS = infP \ ISPPrimeY;
            
            
            Me = wOmega(1) * infPhiIPS(:, 1) * supPhiIPS(:, 1)' * areaOfElement + ...
                wOmega(2) * infPhiIPS(:, 2) * supPhiIPS(:, 2)' * areaOfElement + ...
                wOmega(3) * infPhiIPS(:, 3) * supPhiIPS(:, 3)' * areaOfElement;
            
            Mxe = allp(1, infNodes(1)) * wOmega(1) * infPhiIPS(:, 1) * supPhiIPS(:, 1)' * areaOfElement + ...
                allp(1, infNodes(2)) * wOmega(2) * infPhiIPS(:, 2) * supPhiIPS(:, 2)' * areaOfElement + ...
                allp(1, infNodes(3)) * wOmega(3) * infPhiIPS(:, 3) * supPhiIPS(:, 3)' * areaOfElement;
            
            Mye = allp(2, infNodes(1)) * wOmega(1) * infPhiIPS(:, 1) * supPhiIPS(:, 1)' * areaOfElement + ...
                allp(2, infNodes(2)) * wOmega(2) * infPhiIPS(:, 2) * supPhiIPS(:, 2)' * areaOfElement + ...
                allp(2, infNodes(3)) * wOmega(3) * infPhiIPS(:, 3) * supPhiIPS(:, 3)' * areaOfElement;
            
            
            % Calculate S1, S2, and T
            
            ubarx_ip = ubar(supNodes)*supPhiDxIPS;
            ubary_ip = ubar(supNodes)*supPhiDyIPS;
            ubar_ip = ubar(supNodes)*supPhiIPS;
            
            
            Se1 = wOmega(1) * ubarx_ip(1) * infPhiIPS(:, 1) * supPhiIPS(:, 1)' * areaOfElement + ...
                wOmega(2) * ubarx_ip(2) * infPhiIPS(:, 2) *  supPhiIPS(:, 2)' * areaOfElement + ...
                wOmega(3) * ubarx_ip(3) * infPhiIPS(:, 3) *  supPhiIPS(:, 3)' * areaOfElement;
            
            Se2 = wOmega(1) * ubary_ip(1) * infPhiIPS(:, 1) * supPhiIPS(:, 1)' * areaOfElement + ...
                wOmega(2) * ubary_ip(2) * infPhiIPS(:, 2) *  supPhiIPS(:, 2)' * areaOfElement + ...
                wOmega(3) * ubary_ip(3) * infPhiIPS(:, 3) *  supPhiIPS(:, 3)' * areaOfElement;
            
            Te = wOmega(1) * ubar_ip(1) * infPhiIPS(:, 1) *  supPhiIPS(:, 1)' * areaOfElement + ...
                wOmega(2) * ubar_ip(2) * infPhiIPS(:, 2) *  supPhiIPS(:, 2)' * areaOfElement + ...
                wOmega(3) * ubar_ip(3) * infPhiIPS(:, 3) *  supPhiIPS(:, 3)' * areaOfElement;
            
            
            Ms(infNodes, 1:6) = Ms(infNodes,  1:6) + Me;
            Mxs(infNodes,  1:6) = Mxs(infNodes,  1:6) + Mxe;
            Mys(infNodes,  1:6) = Mys(infNodes,  1:6) + Mye;
            S1s(infNodes,  1:6) = S1s(infNodes,  1:6) + Se1;
            S2s(infNodes,  1:6) = S2s(infNodes,  1:6) + Se2;
            Ts(infNodes,  1:6) = Ts(infNodes,  1:6) + Te;
        end
        
        Msubcell{e} = Ms;
        Mxsubcell{e} = Mxs;
        Mysubcell{e} = Mys;
        Tsubcell{e} = Ts;
        S1subcell{e} = S1s;
        S2subcell{e} = S2s;
        tCell{e} = allt;
        pCell{e} = allp;
    end
    Mcell{part} = Msubcell;
    Mxcell{part} = Mxsubcell;
    Mycell{part} = Mysubcell;
    Tcell{part} = Tsubcell;
    S1cell{part} = S1subcell;
    S2cell{part} = S2subcell;
    subtp{part} = pCell;
    subtt{part} = tCell;
end

toc

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
    Be1 = wOmega(1) * PhiDxIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
        wOmega(2) * PhiDxIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
        wOmega(3) * PhiDxIPS(:, 3) * psiPres(:, 3)' * areaOfElement;
    
    Be2 = wOmega(1) * PhiDyIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
        wOmega(2) * PhiDyIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
        wOmega(3) * PhiDyIPS(:, 3) * psiPres(:, 3)' * areaOfElement;
    
    Be3 = wOmega(1) * PhiIPS(:, 1) *  psiPres(:, 1)' * areaOfElement + ...
        wOmega(2) * PhiIPS(:, 2) *  psiPres(:, 2)' * areaOfElement + ...
        wOmega(3) * PhiIPS(:, 3) * psiPres(:, 3)' * areaOfElement;
    
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
    
    % Viscous terms
    ubarx_ip = ubar(nodes)*PhiDxIPS;
    ubary_ip = ubar(nodes)*PhiDyIPS;
    ubar_ip = ubar(nodes)*PhiIPS;
    
    
    Se1 = wOmega(1) * ubarx_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubarx_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubarx_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
    Se2 = wOmega(1) * ubary_ip(1) * PhiIPS(:, 1) * PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubary_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubary_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
    Te = wOmega(1) * ubar_ip(1) * PhiIPS(:, 1) *  PhiIPS(:, 1)' * areaOfElement + ...
        wOmega(2) * ubar_ip(2) * PhiIPS(:, 2) *  PhiIPS(:, 2)' * areaOfElement + ...
        wOmega(3) * ubar_ip(3) * PhiIPS(:, 3) *  PhiIPS(:, 3)' * areaOfElement;
    
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
UhatNodesp = zeros(numberOfNodes.old, maxWaveNum, numPart);

Btot = @(waveNum) Bs + 1i*waveNum*B3s;
momentumEQNs = @(waveNum) [K + (waveNum.^2)*M - 1i*waveNum*T, sparse(numberOfNodes.new, 2*numberOfNodes.new);
    sparse(numberOfNodes.new,numberOfNodes.new), K + (waveNum.^2)*M - 1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
    S1, S2, K + (waveNum.^2)*M - 1i*waveNum*T];

% Create Preconditioner
Visc = @(waveNum) K + waveNum.^2*M;
Fhat = @(waveNum) blkdiag(-Visc(waveNum), -Visc(waveNum), -Visc(waveNum));

AwnBase = @(waveNum) [-momentumEQNs(waveNum), -Btot(waveNum); -Btot(waveNum)', sparse(numberOfNodes.old,numberOfNodes.old)];
%uStressF = cell(numPart,1);
uStressF = zeros(numPart, 3*numberOfNodes.new, maxWaveNum);
discF = zeros(numPart, 3*numberOfNodes.new, maxWaveNum);
uStressFinf = cell(numPart,1);
%uStressFinf = zeros(numPart, length(p),maxWaveNum);
discFinf = cell(numPart,1);
tic
for i = 1:numPart
    
    
    slicerStress = regStressletF(p',xp(:,i),gammax(i),gammay(i),epsilon, 16*maxWaveNum,L);
    %uStressF{i} = slicerStress(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*10);
    uStressF(i,:,:) = slicerStress(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*16);
    slicerDisc = uDiscF(p', xp(1,i), xp(2,i), gammax(i), gammay(i), 16*maxWaveNum,L);
    discF(i,:,:) = slicerDisc(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*16);
    %discF{i} = slicerDisc(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*10);
    %discF{i}((2*end/3):end,:) = 1i*imag(discF{i}((2*end/3):end,:));
    
    %special points for subtriangle
    allp = [];
    for e = 1:numel(subtp{i})
        allp  = [allp;subtp{i}{e}'];
    end
    [uallp, ~, ic] = unique(allp,'rows');
    
    slicerStressInf = regStressletF(uallp,xp(:,i),gammax(i),gammay(i),epsilon, 16*maxWaveNum,L);
    slicerStressInf = slicerStressInf([ic;ic+length(uallp);ic+2*length(uallp)],:);
    uStressFinf{i} = reshape(slicerStressInf(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*16), E, numel(subtp{i}), 3, maxWaveNum);
    
    
    slicerDiscInf = uDiscF(uallp, xp(1,i), xp(2,i), gammax(i), gammay(i), 16*maxWaveNum,L);
    slicerDiscInf = slicerDiscInf([ic;ic+length(uallp);ic+2*length(uallp)],:);
    discFinf{i} = reshape(slicerDiscInf(:,1:maxWaveNum)/((2*pi)^0.5*maxWaveNum/L*16), E, numel(subtp{i}), 3, maxWaveNum);
end
toc
%compute special int points
RHSopi = @(waveNum,i,e) [-1i*waveNum*Tcell{i}{e}.', sparse(6,2*E);
    sparse(E,6).', -1i*waveNum*Tcell{i}{e}.', sparse(E,6).';
    S1cell{i}{e}.', S2cell{i}{e}.' , -1i*waveNum*Tcell{i}{e}.'];

RHSstressi = @(waveNum,i,e) [+1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.', sparse(6,2*E);
    sparse(E,6).', 1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.', sparse(E,6).';
    sparse(6,2*E) , 1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.'];

RHSstress2xi = @(waveNum,i,e) Re*[1i*waveNum*gammax(i)*Mxcell{i}{e}.', sparse(6,2*E);
    sparse(E,6).', 1i*waveNum*gammax(i)*Mxcell{i}{e}.', sparse(E,6).';
    -gammax(i)*Mcell{i}{e}.', sparse(E,6).', 1i*waveNum*gammax(i)*Mxcell{i}{e}.'];

RHSstress2yi = @(waveNum,i,e) Re*[1i*waveNum*gammay(i)*Mycell{i}{e}.', sparse(6,2*E);
    sparse(E,6).', 1i*waveNum*gammay(i)*Mycell{i}{e}.', sparse(E,6).';
    sparse(E,6).', -gammay(i)*Mcell{i}{e}.', 1i*waveNum*gammay(i)*Mycell{i}{e}.'];

RHSstressWaveScale = @(i,e)[-1i*Tcell{i}{e}.' + 1i*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.' + Re*1i*gammax(i)*Mxcell{i}{e}.'+Re*1i*gammay(i)*Mycell{i}{e}.', sparse(6,2*E);
    sparse(E,6).',   -1i*Tcell{i}{e}.' + 1i*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.' + Re*1i*gammax(i)*Mxcell{i}{e}.' + Re*1i*gammay(i)*Mycell{i}{e}.', sparse(E,6).';
    sparse(6,2*E), -1i*Tcell{i}{e}.' + 1i*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*Mcell{i}{e}.' + Re*1i*gammax(i)*Mxcell{i}{e}.' + Re*1i*gammay(i)*Mycell{i}{e}.'];

RHSstressNoScale = @(i,e) [sparse(6,3*E);
    sparse(6,3*E);
     S1cell{i}{e}.' - Re*gammax(i)*Mcell{i}{e}.',   S2cell{i}{e}.' - Re*gammay(i)*Mcell{i}{e}.', sparse(E,6).'];
 
RHSdisci = @(waveNum, i, e) Re*[-1i*waveNum*(Tcell{i}{e}-us(i)*Mcell{i}{e}).', sparse(6,2*E);
    sparse(E,6).', -1i*waveNum*(Tcell{i}{e}-us(i)*Mcell{i}{e}).', sparse(E,6).';
    S1cell{i}{e}.', S2cell{i}{e}.' , -1i*waveNum*(Tcell{i}{e}-us(i)*Mcell{i}{e}).'];
replacementRHS = cell(numPart,1);

tic
for part = 1:numPart
    replacementRHS{part,maxWaveNum} = zeros(numel(subtp{part}),18);
    
    for e = 1:numel(subtp{part})
        WNmat1 = RHSstressNoScale(part,e);
        WNmat2 = RHSstressWaveScale(part,e);
        for waveIndex = 1:maxWaveNum
            USFi = [uStressFinf{part}(:,e,1,waveIndex); uStressFinf{part}(:,e,2,waveIndex); uStressFinf{part}(:,e,3,waveIndex)];
            DFi = [discFinf{part}(:,e,1,waveIndex); discFinf{part}(:,e,2,waveIndex); discFinf{part}(:,e,3,waveIndex)];
%             
%             replacementRHS{part,waveIndex}(e,:) =  RHSopi(waveNumbers(waveIndex),part,e)*USFi...
%                 + RHSstressi(waveNumbers(waveIndex),part,e)*USFi...
%                 + RHSstress2xi(waveNumbers(waveIndex),part,e)*USFi...
%                 + RHSstress2yi(waveNumbers(waveIndex),part,e)*USFi...
%                 + RHSdisci(waveNumbers(waveIndex),part,e)*DFi;
            %replacementRHS{part,waveIndex}(e,:) =  RHSopi(waveNumbers(waveIndex),part,e).'*USFi;
             replacementRHS{part,waveIndex}(e,:) = (WNmat1*USFi + waveNumbers(waveIndex)*(WNmat2*USFi)) - RHSdisci(waveNumbers(waveIndex),part,e)*DFi;
        end
        %replacementRHS{part,waveIndex} = replacementRHS{part,waveIndex}(ismember(nearbyTris{part},intNearbyTris{part}),:);
        %replacementRHS{part,waveIndex} = reshape(replacementRHS{part,waveIndex},[],1);
    end
end
toc
finalReplacement=cell(numPart,maxWaveNum);
finalNodes = cell(numPart,numPart);
if D>0
    for part = 1:numPart
        finalNodes{part} = unique(t(intNearbyTris{part}(1:end),:));
        for waveIndex = 1:maxWaveNum
            finalReplacement{part,waveIndex} = zeros(numel(finalNodes{part})*3,1);
        end
        for  e = 1:numel(nearbyTris{part})
            nodes = ismember(finalNodes{part},t(nearbyTris{part}(e), :));
            nodes2 = ismember(t(nearbyTris{part}(e), :),finalNodes{part});
            for waveIndex = 1:maxWaveNum
                finalReplacement{part,waveIndex}([nodes;nodes;nodes]) = finalReplacement{part,waveIndex}([nodes;nodes;nodes]) + replacementRHS{part,waveIndex}(e,[nodes2,nodes2,nodes2]).';
            end
        end
    end
end
%%


RHSop = @(waveNum) [-1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new*2);
    sparse(numberOfNodes.new,numberOfNodes.new), -1i*waveNum*T, sparse(numberOfNodes.new,numberOfNodes.new);
    S1, S2 , -1i*waveNum*T];
RHSstress = @(waveNum,i) [+1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*M, sparse(numberOfNodes.new,numberOfNodes.new*2);
    sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*M, sparse(numberOfNodes.new,numberOfNodes.new);
    sparse(numberOfNodes.new, 2*numberOfNodes.new) , 1i*waveNum*(us(i)-Re*gammax(i)*xp(1,i)-Re*gammay(i)*xp(2,i))*M];

RHSstress2x = @(waveNum,i) Re*[1i*waveNum*gammax(i)*Mx, sparse(numberOfNodes.new,numberOfNodes.new*2);
    sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*gammax(i)*Mx, sparse(numberOfNodes.new,numberOfNodes.new);
    -gammax(i)*M, sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*gammax(i)*Mx];

RHSstress2y = @(waveNum,i) Re*[1i*waveNum*gammay(i)*My, sparse(numberOfNodes.new,numberOfNodes.new*2);
    sparse(numberOfNodes.new,numberOfNodes.new), 1i*waveNum*gammay(i)*My, sparse(numberOfNodes.new,numberOfNodes.new);
    sparse(numberOfNodes.new,numberOfNodes.new), -gammay(i)*M, 1i*waveNum*gammay(i)*My];

RHSdisc = @(waveNum, i) Re*[-1i*waveNum*(T-us(i)*M), sparse(numberOfNodes.new,numberOfNodes.new*2);
    sparse(numberOfNodes.new,numberOfNodes.new), -1i*waveNum*(T-us(i)*M), sparse(numberOfNodes.new,numberOfNodes.new);
    S1, S2 , -1i*waveNum*(T-us(i)*M)];

% discF = reshape(cell2mat(discF), 3*numberOfNodes.new, numPart, maxWaveNum);
% uStressF = reshape(cell2mat(uStressF),  3*numberOfNodes.new, numPart, maxWaveNum);
discF = permute(discF, [2,1,3]);
uStressF = permute(uStressF, [2,1,3]);
%rank one update to increase rank
u = [zeros(3*numberOfNodes.new,1); ones(numberOfNodes.old,1)];
u = sparse(u);

% Dirichlet boundary
Dirichlet123 = [Dirichlet; Dirichlet+numberOfNodes.new;Dirichlet+2*numberOfNodes.new];
parfor waveIndex = 2:maxWaveNum
    
    
    
    [FhatL, FhatU, FhatP, FhatQ] = lu(Fhat(waveNumbers(waveIndex)));
    schurCompi = 1*mean(mean(abs(Bs)))*eq(waveIndex,1)*ones(numberOfNodes.old)+Mp;
    [schurL, schurU] = lu(schurCompi);
    tol = 10^(-7);
    maxit = 50;
    buttshit = zeros(6, 3, numPart); %you need to create this because of how matlab does slicing
    %buttshit2 = zeros(numberOfNodes.new, numPart); %you need to create this because of how matlab does slicing
    RHSbase = RHSop(waveNumbers(waveIndex));
    F = zeros(wnl,1);
    ABase = AwnBase(waveNumbers(waveIndex));
    A=ABase;
    A(Dirichlet123, :) = 0;
    uStressFWave = uStressF(:,:,waveIndex);
    discFWave = discF(:,:,waveIndex);
    for i = 1:numPart
        A(1:a1/3,1:a1/3) = ABase(1:a1/3,1:a1/3)-1i*waveNumbers(waveIndex)*us(i)*M;
        A((a1/3+1):(2*a1/3),(a1/3+1):(2*a1/3)) = A(1:a1/3,1:a1/3);
        A((2*a1/3+1):(a1),(2*a1/3+1):(a1)) = A(1:a1/3,1:a1/3);
        
        F(:,i) = [RHSbase*uStressFWave(:,i)...
            + RHSstress(waveNumbers(waveIndex),i)*uStressFWave(:,i)...
            + RHSstress2x(waveNumbers(waveIndex),i)*uStressFWave(:,i)...
            + RHSstress2y(waveNumbers(waveIndex),i)*uStressFWave(:,i)...
            - RHSdisc(waveNumbers(waveIndex),i)*discFWave(:,i);
            sparse(numberOfNodes.old,1)];
        if D>0
            replaceNodes = [finalNodes{i},finalNodes{i} + numberOfNodes.new,finalNodes{i} + 2*numberOfNodes.new];
            F(replaceNodes,i) = finalReplacement{i,waveIndex};
        end
        F(Dirichlet123, i) = (-uStressFWave(Dirichlet123,i) + Re*discFWave(Dirichlet123,i));
        
        if waveIndex == 1
            A = A + 100*mean(mean(abs(Bs)))*(u*u');
        end
        
        Uhat = gmres(A, F(:, i), 25, tol,maxit, @(x) PCbackSolve(x,A(1:a1, (a1+1):end),FhatL,FhatU, FhatP, FhatQ, schurL, schurU, a1));
        
        buttshit(:,1, i) = Uhat(particleNodes{i});
        buttshit(:,2, i) = Uhat(particleNodes{i} + numberOfNodes.new);
        %buttshit(:,3, i) = Uhat(particleNodes{i} + 2*numberOfNodes.new);
        %buttshitp(:, i) = Uhat((1:numberOfNodes.old) + 3*numberOfNodes.new);
    end
    UhatNodes1(:, waveIndex, :) = buttshit(:, 1, :);
    UhatNodes2(:, waveIndex, :) = buttshit(:, 2, :);
    %UhatNodes3(:, waveIndex, :) = buttshit(:, 3, :);
    %UhatNodesp(:, waveIndex, :) = buttshitp(:, :);
    
end

%inverse fourier transform and evaluate velocity

for i = 1:numPart
    %UhatNodes2(1:6,:,i) = UhatNodes2(particleNodes{i},:,i);
    %UhatNodes2(1:6,1,i) = UhatNodes2(:,1,i);
    %UhatNodes1(1:6,:,i) = UhatNodes1(particleNodes{i},:,i);
    %UhatNodes1(1:6,1,i) = UhatNodes1(:,1,i);
end
UhatNodes2 = UhatNodes2(1:6,:,:);
UhatNodes1 = UhatNodes1(1:6,:,:);
 Unodes2 = sum(real(UhatNodes2), 2)/maxWaveNum;
 Unodes1 = sum(real(UhatNodes1), 2)/maxWaveNum;
%Unodes2 = real(UhatNodes2)/maxWaveNum;
%Unodes1 = real(UhatNodes1)/maxWaveNum;

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
    
    velocity(i,1) = Unodes1(:,i)' * PhiIPS/Re;
    velocity(i,2) = Unodes2(:,i)' * PhiIPS/Re;
end