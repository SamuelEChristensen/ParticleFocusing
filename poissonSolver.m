function [U,pOld,tOld]=poissonSolver(p,t,f,fb)

%keeping fb because it will be a pain to put back in. it should always be
%zero

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
Dirichlet = [];
for i=1 : length(Dirichlet_e)
    nodes = Dirichlet_e(i, :);
    Dirichlet = [Dirichlet; [nodes(1), S(nodes(1), nodes(2)), nodes(2)]];
end
Dirichlet=unique(Dirichlet);
Fb=zeros(length(Dirichlet),1);
Fb(1:length(Dirichlet))=reshape(fb(p(:,Dirichlet)'),length(Dirichlet),1);


% Initialisation of K, F
K = sparse(numberOfNodes.new,numberOfNodes.new);
F = zeros(numberOfNodes.new,1);
 
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
     
    
    K(nodes, nodes) = K(nodes, nodes) + Ke;
     
       
    Fe = wOmega(1) * PhiIPS(:, 1) * f(ip(1, :)) * areaOfElement + ...
         wOmega(2) * PhiIPS(:, 2) * f(ip(2, :)) * areaOfElement + ...
         wOmega(3) * PhiIPS(:, 3) * f(ip(3, :)) * areaOfElement; 
     
    F(nodes) = F(nodes) + reshape(Fe, length(nodes), 1);   

end


% Dirichlet boundary

K(Dirichlet, :) = 0;
K(Dirichlet, Dirichlet) = -speye(numel(Dirichlet));

F(Dirichlet) = Fb(1:numel(Dirichlet));



U = -(K\F);


tol=10^(-10);
maxit=100;
%U=gmres(A,F,30,tol,maxit, PCG);