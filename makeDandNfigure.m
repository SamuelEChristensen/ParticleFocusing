X = load('lotsOfPoints');
%this does not work ATM
paramSetCts = {[0.02  0.000825  150 5]};
[a,b] = meshgrid(1:4,1:4);
paramSetNandD = [a(:),b(:)];
xp = [2; 0.275];
velocities = zeros(3,length(paramSetSing), size(xp,2),2);
numTris = zeros(3,length(paramSetSing));

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
figure
 parami = paramSetCts{1};
 fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
for i = 1:length(paramSetNandD)
   

numTris(3,i) = length(t);
maxWaveNum = parami(3);
L = parami(4);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;



[bgFlow,~,~] = poissonSolver(p,t,f,fb);

[ucts,~,~] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L,paramSetNandD(i,1), paramSetNandD(i,2));
velocities(3,i,:,:) = (2*pi)^0.5*maxWaveNum/L*ucts*2;
end
fucky = reshape(velocities(3,:,2),size(a));
fuckx = reshape(velocities(3,:,1),size(a));
ass = ((fucky-fucky(4,4)).^2+fuckx.^2).^0.5/fucky(4,4);