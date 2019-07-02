%test case for oseenSolver with gaussian approximation of delta function
fd=@(p) sqrt(sum(p.^2,2))-1;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
[p,t]=distmesh2d(fd,@huniform,0.04,[-1,-1;1,1],[0  0]);
maxWaveNum = 128;
L = 10;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
bestGuess = zeros(1,6);
parfor i = 1:6
    xp = [0.15*i; 0];
    [u,U,pold,told]=velocitySolve(p,t, 0.05, waveNumbers, xp);
    bestGuess(i) = u(1);
    
end
%figure
%plot(bestGuess)