%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = {[0.01   0.01  512]  [0.0141   0.01  256]   [0.028   0.01  512]    [0.0141   0.01  512]};
paramSet = {[0.02   0.0001  200]};
%paramSet = {[0.005   0.002  1500] };
xp = [2*ones(size(0.1:0.1:0.45)); 0.1:0.1:0.45];
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
L = 4;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow');
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end

 X = readmatrix('hinchComsol');
 Y = readmatrix('hinchData');


figure
hold on
plot(xp(2,:), real(velocities(:,:,2)))
plot(X(:,1),X(:,2))
plot(Y(:,1),Y(:,2))