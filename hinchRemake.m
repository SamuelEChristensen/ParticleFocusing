%test case for oseenSolver with gaussian approximation of delta function

fd=@(p) drectangle(p,0,4,0,1);
[p,t]=distmesh2d(fd,@huniform,0.045,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
maxWaveNum = 64*3;
L = 6;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
%moviePlot = cell(30, 1);
tic

xp = [2*ones(size(0.1:0.05:0.45)); 0.1:0.05:0.45];
[u,pold,told]=velocitySolveExperimental(p,t, 0.035, waveNumbers, xp);
%moviePlot{i} = U;
velocities = (2*pi)^0.5*maxWaveNum/L*u;
(2*pi)^0.5*maxWaveNum/L*u
toc

X = readmatrix('hinchComsol');
Y = readmatrix('hinchData');


figure
hold on
plot(xp(2,:), real(velocities(:,2)))
plot(X(:,1),X(:,2))
plot(Y(:,1),Y(:,2))
