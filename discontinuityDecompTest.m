%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
paramSet = {[0.025  0.0015  200 5  2] };
%paramSet = {[0.005   0.002  1500]};
xp = [2*ones(size(0.1:0.05:0.45)); 0.1:0.05:0.45];
%xp = [2;0.2];
velocities = zeros(length(paramSet), size(xp,2),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));
i=1;
parami = paramSet{i};
fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
%fh=@(p) parami(1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
tic
fb=@(p) sol(p);
figure
for i = 1:length(paramSet)
    parami = paramSet{i};
    toc
    maxWaveNum = parami(3);
    L = parami(4);
    
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    tic
    sol=@(x) zeros(size(x(:,1)));
    f=@(x) -ones(size(x(:,1)));
    
    
    fb=@(p) sol(p);
    [bgFlow,~,~] = poissonSolver(p,t,f,fb);
    [u,pold,told] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L,1,2,1);
    velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
    toc
end

X = xlsread('hinchComsol_75');
Y = xlsread('hinchData_75');
figure
hold on
plot(X(:,1),X(:,6));
plot(Y([1,2,4:9],1),Y([1,2,4:9],2));
plot(xp(2,:), real(velocities(:,:,2)));
names = {'Comsol \epsilon=10^{-3}','Hinch'};
for i = 1:numel(paramSet)
    names(2+i) = {['inFocus Param set ',num2str(i)]};
end
title('1-D Poiseuille migration velocity')
legend(names)