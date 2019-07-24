%test case for oseenSolver with gaussian approximation of delta function

%notes for future improvement. Calculate all solves for one wave number at
%once, then move on to next wave number. This may help with
%parrelelization.
fd=@(p) sqrt(sum(p.^2,2))-1;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
[p,t]=distmesh2d(fd,@huniform,0.05,[-1,-1;1,1],[0  0]);
maxWaveNum = 512;
L = 5;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
%moviePlot = cell(30, 1);
tic

xp = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8; 0,0,0,0,0,0,0,0];
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
plot(xp(1,:), 16*real(velocities(:,1)))
plot(ass(:,1),ass(:,11))
