X = load('lotsOfPoints');

paramSetSing = {[0.1  0.05  400 5]  [0.05  0.025  400  5]  [0.025  0.0125 400  5]   [0.0125  0.0065 400  5]   [0.0075  0.004 600  5]   [0.005  0.0025 800  5]};
paramSetDisc = {[0.1  0.025 300 5]   [0.05  0.025  300 5]   [0.025  0.0125  300 5]    [0.0125  0.00625  300 5]    [0.0075  0.004  400 5]   [0.005  0.0025 400  5]};
paramSetCts = {[0.1  0.005  150 5]   [0.05  0.0025  150 5]   [0.025  0.00125  150 5]    [0.0125  0.000625  150 5]    [0.0075  0.0005  150 5]    [0.005  0.0003 150  5]};
xp = [2; 0.3];
velocities = zeros(3,length(paramSetSing), size(xp,2),2);
numTris = zeros(3,length(paramSetSing));


gridSize = [0.1 0.05 0.025  0.0125  0.0075, 0.005];
sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
figure
for i = 1:length(paramSetSing)
    parami = paramSetSing{i};
 fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
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
 fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
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
 fd=@(p) drectangle(p,0,4,0,1);
fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
[p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
numTris(3,i) = length(t);
maxWaveNum = parami(3);
L = parami(4);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;



[bgFlow,~,~] = poissonSolver(p,t,f,fb);

[ucts,~,~] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L, 2, 2);
velocities(3,i,:,:) = (2*pi)^0.5*maxWaveNum/L*ucts*2;
end


%Y = xlsread('hinchData');

errorSing = zeros(length(paramSetSing),1);
errorDisc = zeros(length(paramSetDisc),1);
errorCts = zeros(length(paramSetCts),1);
 for i = 1:(length(paramSetSing)-1)
errorSing(i) = mean(abs(velocities(1,end,:,2)-velocities(1,i,:,2)))./abs(velocities(3,end,:,2));
 end
 for i =1:(length(paramSetDisc)-1)
errorDisc(i) =  mean(abs(velocities(2,end,:,2)-velocities(2,i,:,2)))./abs(velocities(3,end,:,2));
 end
 for i=1:(length(paramSetCts)-1)
errorCts(i) =  mean(abs(velocities(3,end,:,2)-velocities(3,i,:,2)))./abs(velocities(3,end,:,2));
 end

figure
loglog(gridSize, errorSing,'-*', gridSize, errorDisc,'-square',gridSize, errorCts,'-o','lineWidth',3.0)
%title('Self Convergence of 3 Different Formulations')
XTickLabels = flipud(cellstr(num2str(round(log2(gridSize(1:(end-1)))'), '2^{%d}')));
YTickLabels = (cellstr(num2str([-11,-8,-5,-2]', '2^{%d}')));
ax = gca;
ax.FontSize = 13; 
set(gca,'xticklabel',XTickLabels,'xtickmode','manual');
set(gca,'yticklabel',YTickLabels,'ytickmode','manual');
xticks(fliplr([2.^round(log2(gridSize(1:(end-1))))]))
yticks([2^(-11),2^(-8),2^(-5),2^(-2)])
xlabel('Grid Size','fontSize',22)
ylabel('Relative Error','fontSize',22)
% 
% plot(t,Y_init_fit,'color',color,'LineWidth',3);
% color = hsv2rgb([0 0.5 1]);
% plot(t,Y_1,'color',color,'LineWidth',3);
% color = hsv2rgb([0 0.75 1]);
% plot(t,Y_2,'color',color,'LineWidth',3);