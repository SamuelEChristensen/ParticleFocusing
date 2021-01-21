%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = {[0.01   0.01  512]  [0.0141   0.01  256]   [0.028   0.01  512]    [0.0141   0.01  512]};
paramSet = {[0.01  0.025  160  4] };
paramSet = {[0.1  0.05  300 5]   [0.05  0.025  300 5]   [0.025  0.0125  300 5]    [0.0125  0.00625  300 5]    [0.0075  0.0025  400 5]};
paramSet = {[0.01   0.0025  400  4] };
xp = [2*ones(size(0.1:0.05:0.45)); 0.1:0.05:0.45];
%xp = [2;0.2];
velocities = zeros(length(paramSet), size(xp,2),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for i = 1:length(paramSet)
    parami = paramSet{i};
      
fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(parami(1)+1.4*drectangle(p, 1.99, 2.01, 0.1, 0.5).^3,0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
maxWaveNum = parami(3);
L = parami(4);
% z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
% z = circshift(z, -maxWaveNum/2+1);
%waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told] = velocitySolveStressletDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end


X = xlsread('hinchComsol_75');
Y = xlsread('hinchData_75');
figure
hold on
% plot(X(:,1),X(:,2));
% plot(Y([1,2,4:9],1),Y([1,2,4:9],2));
 plot(X(:,1),X(:,6))
 plot(Y(:,1),Y(:,2))
plot(xp(2,:), real(velocities(:,:,2))/75)