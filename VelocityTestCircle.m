%test case for oseenSolver with gaussian approximation of delta function
fd=@(p) sqrt(sum(p.^2,2))-1;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
[p,t]=distmesh2d(fd,@huniform,0.075,[-1,-1;1,1],[0  0]);
maxWaveNum = 256;
L = 10;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
velocities = zeros(1, 30);
%moviePlot = cell(30, 1);
tic
for i = 5
    xp = [0.5; 0];
    [u,U,pold,told]=velocitySolve(p,t, 0.05, waveNumbers, xp);
    %moviePlot{i} = U;
    velocities(i) = (2*pi)^0.5*maxWaveNum/L*u(1);
    (2*pi)^0.5*maxWaveNum/L*u
    toc
end
%figure
%plot(bestGuess)