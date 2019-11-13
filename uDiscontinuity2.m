%reccomend doing on computer with lots of memory
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = { [0.0075   0.005  180]};
paramSet = {[0.005   0.0025  220]};
xp = [2; 0.25];%[2*ones(size(0.1:0.05:0.45)); 0.1:0.05:0.45];
velocities = cell(length(paramSet), length(xp),2);
Ux = zeros(1000,length(paramSet));
Uy = zeros(1000,length(paramSet));
for i = 1:length(paramSet)
    parami = paramSet{i};
    
fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1), parami(1)+1.4*drectangle(p, 1.99, 2.01, 0.1, 0.5).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
maxWaveNum = parami(3);
L = 4;
% z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
% z = circshift(z, -maxWaveNum/2+1);
%waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[u,p,t]=velocitySolveExperimental(p,t, parami(2), waveNumbers, xp);
    for j = 1:length(1)
        velocities(i,j,1) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,1)*2};
        velocities(i,j,2) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,2)*2};
    end   
toc
    for j = 1:1000
        xpi = [2; 0.5*j/1000];
        xpn = whatTriangleIsThisPointIn(p,t,xpi);
        nodes = t(xpn, :);
        %particleNodes{j} = nodes;
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
 
    
    %bgFlowIPS = bgFlow(nodes);
    for k = 1:length(1)
        Ux(j,i) = velocities{i,k,1}(nodes).' * PhiIPS;
        Uy(j,i) = velocities{i,k,2}(nodes).' * PhiIPS;
    end
    %gammax(j) = bgFlowIPS * PhiDxIPS;
    %gammay(j) = bgFlowIPS * PhiDyIPS;
    xvec = 0.5/1000:0.5/1000:0.5;
    end
    figure
    for j=1:length(paramSet)
    subplot(2,length(paramSet),j)
    plot(xvec,Uy(:,j)-5/9*(8*(xp(2,1)-0.5))^2*heaviside(xvec'-xp(2,1)))
    subplot(2,length(paramSet),j+length(paramSet))
    plot(xvec,Uy(:,j))
    end
end
