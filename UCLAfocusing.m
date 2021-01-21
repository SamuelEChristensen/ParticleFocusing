
paramSet = {[0.014   0.0015  150]};
%  smallest grid size epsilon  max wavenumber (just keep it at like 300)
 

% this for cross section plot
[X,Y] =  meshgrid(-0.9:0.1:5);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
xp=xp(:,fd(xp')<-0.05); %mess here
%velocities = zeros(length(paramSet), length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for j = 5:8
xpj = xp(:, (1:56)+(j-1)*56);
for i = 1:length(paramSet)
    parami = paramSet{i};
  

maxWaveNum = parami(3);
L = 4;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told]=velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xpj, bgFlow',4,2);
velocities(i,(1:56)+(j-1)*56,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
end
% 
figure
hold on
polyU = polyshape(logoVertsU(orderU,1),logoVertsU(orderU,2));
polyC = polyshape(logoVertsC(orderC,1),logoVertsC(orderC,2));
polyL = polyshape(fuck(:,1),fuck(:,2));
polyA = polyshape(logoVertsA(orderA,1),logoVertsA(orderA,2));
plot(polyU)
plot(polyC)
plot(polyL)
plot(polyA)
v = [-1 -1; -1 0.5; 4 0.5; 4  -1];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none');

quiver(xp(1,:),xp(2,:),velocities(1,:,1),velocities(1,:,2),'b')
%quiver(xp(1,:),xp(1,:),-velocities(1,:,2),velocities(1,:,1),'b')

figure
hold on
quiver(xp(1,:),xp(2,:),velocities(1,:,1),velocities(1,:,2),'b')
%[p,t]=distmesh2d(fd,@huniform,0.03 ,bbox,pfix);
Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', velocities(1,:,2)','natural');
Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', velocities(1,:,1)','natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.09*max(xlength,ylength);
sp=p(ix,:);
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,600], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2))
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.2,'ByRows',1,'DataScale',1);
ix2=fd(eqPoints)<0;
figure
hold on
plotFESol(p,t,(Fx(p).^2+Fy(p).^2).^0.25/1.5)
scatter3(eqPoints(ix2,1),eqPoints(ix2,2),repmat(0.5,length(eqPoints(ix2,1)),1),50,'k','filled')
axis image
set(gca,'visible','off')