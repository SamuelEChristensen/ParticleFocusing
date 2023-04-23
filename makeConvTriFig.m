paramSetSing = {[0.1  0.05  256 5]  [0.05  0.025  256  5]  [0.025  0.0125 256  5]   [0.0125  0.0065 256  5]   [0.006125  0.004 256  5]  [0.0035  0.002 256  5]};
paramSetDisc = {[0.1  0.025 256 5]   [0.05  0.0125  256 5]   [0.025  0.00625  256 5]    [0.0125  0.00325  256 5]    [0.006125  0.0016  256 5] [0.0035  0.0008 256  5]};
paramSetCts = {[0.1  0.005  128 5]   [0.05  0.00025  128 5]   [0.025  0.000125  128 5]    [0.0125  0.000625  128 5]    [0.006125  0.0005  128 5]    [0.0035  0.0003 128  5]};
xp = [-0.2; -0.1];
velocities = zeros(3,length(paramSetSing), size(xp,2),2);
numTris = zeros(3,length(paramSetSing));
pv = [-0.51  -0.51; -0.51  0.51; 0.51  0.51;  0.51  -0.51;  -0.51  -0.51];

fd=@(p) drectangle(p,-0.51,0.51,-0.51,0.51);

gridSize = [0.1 0.05 0.025  0.0125  0.006125 0.0035];
sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
figure
for i = 1:length(paramSetSing)
    parami = paramSetSing{i};
    fh=@(p) min(parami(1)+max(0,0.5*drectangle(p, -0.19, 0.-0.21, -0.09, -0.11).^3),0.1);  % custom mesh distance that puts more points in one quadrant
    
    [p,t]=distmesh2d(fd,fh,parami(1),[-0.51, -0.51;0.51, 0.51],pv);
    numTris(1,i) = length(t);
    maxWaveNum = parami(3);
    L = parami(4);
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    
    
    
    [bgFlow,~,~] = poissonSolver(p,t,f,fb);
    
    [using,pold,told]=velocitySolveFull(p,t, parami(2), waveNumbers, xp, bgFlow');
    velocities(1,i,:,:) = (2*pi)^0.5*maxWaveNum/L*using*2;
    
     
    parami = paramSetDisc{i};
    maxWaveNum = parami(3);
    L = parami(4);
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    [udisc,~,~] = velocitySolveStressletDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L);
    velocities(2,i,:,:) = (2*pi)^0.5*maxWaveNum/L*udisc*2;
    
    
    parami = paramSetCts{i};
    maxWaveNum = parami(3);
    L = parami(4);
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    [ucts,~,~] = velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L,1,2, 50);
    velocities(3,i,:,:) = (2*pi)^0.5*maxWaveNum/L*ucts*2;
end



errorSing = zeros(length(paramSetSing),1);
errorDisc = zeros(length(paramSetDisc),1);
errorCts = zeros(length(paramSetCts),1);
for i = 1:(length(paramSetSing)-1)
    errorSing(i) = mean(abs(velocities(1,end,:,1)-velocities(1,i,:,1)))./abs(velocities(3,end,:,1));
end
for i =1:(length(paramSetDisc)-1)
    errorDisc(i) =  mean(abs(velocities(2,end,:,1)-velocities(2,i,:,1)))./abs(velocities(3,end,:,1));
end
for i=1:(length(paramSetCts)-1)
    errorCts(i) =  mean(abs(velocities(3,end,:,1)-velocities(3,i,:,1)))./abs(velocities(3,end,:,1));
end

figure
loglog(gridSize, errorSing,'-^', gridSize, errorDisc,'-square',gridSize, errorCts,'-o','lineWidth',3.0)
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