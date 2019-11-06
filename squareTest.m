%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = {[0.01   0.01  512]  [0.0141   0.01  256]   [0.028   0.01  512]    [0.0141   0.01  512]};
paramSet = {[0.014   0.01  200]   [0.02  0.014  150]};
%paramSet = {[0.005   0.002  1500] };
[X,Y] =  meshgrid(-.4:0.05:0);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
velocities = zeros(length(paramSet), length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for i = 1:length(paramSet)
    parami = paramSet{i};
  
fd=@(p) drectangle(p,-0.5,0.5,-0.5,0.5);
fh=@(p) min(parami(1)+max(0,0.5*drectangle(p, -0.5, 0, -0.5, 0).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[-0.5, -0.5;0.5, 0.5],[-0.5  -0.5; 0.5  -0.5;  0.5  0.5;  -0.5   0.5]);
maxWaveNum = parami(3);
L = 4;
% z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
% z = circshift(z, -maxWaveNum/2+1);
%waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told]=velocitySolveFull(p,t, parami(2), waveNumbers, xp, bgFlow');
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
% 
 figure
quiver(X',Y',velocities(1,:,1),velocities(1,:,2))