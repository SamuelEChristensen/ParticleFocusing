
paramSet = {[0.014   0.001  150]};
%  smallest grid size epsilon  max wavenumber (just keep it at like 150)
 

% this for cross section plot
%[X,Y] =  meshgrid(-0.9:0.05:5);
[X,Y] =  meshgrid(-0.9:0.04:5,-0.9:0.07:5);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
xp=xp(:,fd(xp')<-0.01); %mess here
[X,Y] =  meshgrid(1:0.02:1.2,0.5:0.035:0.8);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp2 = [X,Y]';
xp2=xp2(:,fd(xp2')<-0.005); %mess here
xp = [xp,xp2];
velocities = zeros(length(paramSet), length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for i = 1:length(paramSet)
    parami = paramSet{i};
  

maxWaveNum = parami(3);
L = 4;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told]=velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',4,1,2);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
% 
fullXP = xp';
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = fuck;
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.02;
sp=p(ix,:);
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,300], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2))
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
finalIC = zeros(length(p),1);
finalIC(ix) = ic;
figure
plotFESol(p,t,finalIC)
hold on
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),140,'k','filled')
plot(wah(:,1)/1787,wah(:,2)/1787)
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)


%figure
%hold on
for i = 1:5:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,100], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2),'b')
arowscal = 1;
quiver(Y(:,1), Y(:,2), Fx(Y)/5, Fy(Y)/5,'b','LineWidth',1.5)
end
% [eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
% 
% hold on
% scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),140,'k','filled')
% plot(wah(:,1)/1787,wah(:,2)/1787)
% set(gca,'xtick',[],'ytick',[])
% xlabel('x','fontSize',24)
% ylabel('y','fontSize',24)