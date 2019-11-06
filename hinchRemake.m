%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
paramSet = {[0.03   0.04  150]  [0.03   0.03  150]   [0.02   0.02  150]  [0.01   0.015  180] [0.01   0.01  180]    [0.0075   0.005  180]};
%paramSet = {[0.01   0.005  178]};
%paramSet = {[0.005   0.002  1500] };
xp = [2*ones(size(0.1:0.05:0.45)); 0.1:0.05:0.45];
velocities = zeros(length(paramSet), length(xp),2);
for i = 1:length(paramSet)
    parami = paramSet{i};
    
fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(parami(1)+1.4*drectangle(p, 1.99, 2.01, 0.1, 0.5).^3,0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
maxWaveNum = parami(3);
L = 4;
% z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
% z = circshift(z, -maxWaveNum/2+1);
%waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[u,pold,told]=velocitySolveExperimental(p,t, parami(2), waveNumbers, xp);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
% 
 X = readmatrix('hinchComsol');
 Y = readmatrix('hinchData');


figure
hold on
plot(xp(2,:), real(velocities(:,:,2)))
plot(X(:,1),X(:,2))
plot(Y(:,1),Y(:,2))

% 
% X = readmatrix('hinchComsol_75');
% Y = readmatrix('hinchData_75');
% figure
% hold on
% plot(xp(2,:), real(velocities(:,:,2))/75)
% plot(X(:,1),X(:,6))
%  plot(Y(:,1),Y(:,2))
