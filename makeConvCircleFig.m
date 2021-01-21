X = load('lotsOfPoints');

paramSetSing = {[0.1  0.05  400 5]  [0.05  0.025  400  5]  [0.025  0.0125 400  5]   [0.0125  0.0065 400  5]   [0.009  0.005 400  5]};
paramSetDisc = {[0.1  0.025 300 5]   [0.05  0.0125  300 5]   [0.025  0.00625  300 5]    [0.0125  0.003125  300 5]  [0.009  0.0025  300 5]};
paramSetCts = {[0.1  0.005  200 5]   [0.05  0.0025  200 5]   [0.025  0.00125  200 5]    [0.0125  0.000625  200 5]   [0.009  0.000425  200 5]};
xp = [0; 0.15];
velocities = zeros(3,length(paramSetSing), size(xp,2),2);
numTris = zeros(3,length(paramSetSing));
gridSize = [0.1  0.05  0.025  0.0125  0.009 ];
sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
figure
for i = 1:length(paramSetSing)
    parami = paramSetSing{i};
fd=@(p) sqrt(sum(p.^2,2))-0.7;
[p,t]=distmesh2d(fd,@huniform,parami(1), [-1, -1; 1, 1],[]);
numTris(1,i) = length(t);
maxWaveNum = parami(3);
L = parami(4);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;



[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[using,pold,told]=velocitySolveFull(p,t, parami(2), waveNumbers, xp, bgFlow');

velocities(1,i,:,:) = (2*pi)^0.5*maxWaveNum/L*using*2;
end
for i = 1:length(paramSetDisc)
    parami = paramSetDisc{i};
fd=@(p) sqrt(sum(p.^2,2))-0.7;
[p,t]=distmesh2d(fd,@huniform,parami(1), [-1, -1; 1, 1],[]);
numTris(2,i) = length(t);
maxWaveNum = parami(3);
L = parami(4);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;



[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[udisc,~,~] = velocitySolveStressletDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L);
velocities(2,i,:,:) = (2*pi)^0.5*maxWaveNum/L*udisc*2;
end

for i = 1:length(paramSetCts)
    parami = paramSetCts{i};
fd=@(p) sqrt(sum(p.^2,2))-0.7;
[p,t]=distmesh2d(fd,@huniform,parami(1), [-1, -1; 1, 1],[]);
numTris(3,i) = length(t);
maxWaveNum = parami(3);
L = parami(4);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;



[bgFlow,~,~] = poissonSolver(p,t,f,fb);

[ucts,~,~] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L, 1, 2);
velocities(3,i,:,:) = (2*pi)^0.5*maxWaveNum/L*ucts*2;
end


%Y = xlsread('hinchData');

errorSing = zeros(length(paramSetSing),1);
errorDisc = zeros(length(paramSetDisc),1);
errorCts = zeros(length(paramSetCts),1);
 for i = 1:(length(paramSetSing)-1)
errorSing(i) = mean(abs(velocities(1,end,:,2)-velocities(1,i,:,2)))./abs(velocities(1,end,:,2));
 end
 for i =1:(length(paramSetDisc)-1)
errorDisc(i) =  mean(abs(velocities(2,end,:,2)-velocities(2,i,:,2)))./abs(velocities(2,end,:,2));
 end
 for i=1:(length(paramSetCts)-1)
errorCts(i) =  mean(abs(velocities(3,end,:,2)-velocities(3,i,:,2)))./abs(velocities(3,end,:,2));
 end

figure
loglog(gridSize, errorSing,'-*', gridSize, errorDisc,'-square',gridSize, errorCts,'-o')
title('Self Convergence of 3 Different Formulations')
xlabel('# of triangles')
ylabel('Relative Error')
legend('Singular','Discontinuous','Continuous')