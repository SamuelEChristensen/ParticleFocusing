%test case for oseenSolver with gaussian approximation of delta function

%fd=@(p) sqrt(sum(p.^2,2))-1;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
%[p,t]=distmesh2d(fd,@huniform,0.03,[-1,-1;1,1],[0  0]);
% fh=@(p) min(0.001+1*abs(dcircle(p,0.5,0,0)).^2,0.1);
% [p,t]=distmesh2d(fd,fh,0.001,[-1,-1;1,1],[0,0]);
fh=@(p) min(0.04+2*drectangle(p,-0.01,0.01,0,1).^3,0.1);
fd=@(p) sqrt(sum(p.^2,2))-2;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
[p,t]=distmesh2d(fd,fh,0.04,[-1,-1;1,1],[0  0]);
maxWaveNum = 256*2;
L = 4;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
%moviePlot = cell(30, 1);
tic

xp = [ zeros(size(0.1:0.1:0.9));0.1:0.1:0.9];
 profile on
[u,pold,told]=velocitySolveExperimental(p,t, 0.03, waveNumbers, xp);
 profile viewer
 profile off
%moviePlot{i} = U;
velocities = (2*pi)^0.5*maxWaveNum/L*u;
(2*pi)^0.5*maxWaveNum/L*u
toc

X = readmatrix('circleComsol');

fuck = (X(:,7) ==0);
ass = X(fuck,:);
figure
hold on
plot(xp(2,:), real(velocities(:,2)))
plot(ass(:,1),ass(:,11))
