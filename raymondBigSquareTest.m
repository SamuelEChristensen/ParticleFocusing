%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = {[0.01   0.01  512]  [0.0141   0.01  256]   [0.028   0.01  512]    [0.0141   0.01  512]};
paramSet = {[0.1   0.05  175]};
%  smallest grid size epsilon  max wavenumber (just keep it at like 300)
 
fd = @(p) drectangle(p,-5,5,-5,5);
fh = @(p) min(parami(1)+max(0,0.5*drectangle(p, 0, 5, 0, 5).^3),0.5);
bgFunc = @(x,y) x+x.^2;
% this for cross section plot
% [X,Y] =  meshgrid(0:1:5);
% X = reshape(X,numel(X),1);
% Y = reshape(Y,numel(Y),1);
% xp = [X,Y]';
% xp=xp(:,fd(xp')<-0.03); %mess here
xp=[0;0];
velocities = zeros(length(paramSet), 1,2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for i = 1:length(paramSet)
    parami = paramSet{i};
  
fh = @(p) min(parami(1)+max(0,0.5*drectangle(p, 0, 5, 0, 5)),0.5);
[p,t]=distmesh2d(fd,fh,parami(1),[-5,-5;5,5],[5,-5;5,5;-5,5;-5,-5]);
maxWaveNum = parami(3);
L = 4;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
Re=1;
[u,pold,told]=velocitySolveRaymond(p,t, parami(2), waveNumbers, xp, bgFunc,Re);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
% 
%figure
%x = 0:0.01:5;
%squareaxes=gca;

%axes(squareaxes)
%hold on
%quiver(xp(2,:),xp(1,:),velocities(1,:,2),velocities(1,:,1),'b')

