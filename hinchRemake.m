%test case for oseenSolver with gaussian approximation of delta function

fd=@(p) drectangle(p,0,4,0,1);
[p,t]=distmesh2d(fd,@huniform,0.06,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
maxWaveNum = 256;
L = 3;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
%moviePlot = cell(30, 1);
tic

xp = [2*ones(size(0.1:0.025:0.4)); 0.1:0.025:0.4];
 profile on
[u,pold,told]=velocitySolve(p,t, 0.05, waveNumbers, xp);
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
plot(xp(2,:), 1/(4*0.11^3)*real(velocities(:,2)))
plot(ass(:,1),ass(:,11))
