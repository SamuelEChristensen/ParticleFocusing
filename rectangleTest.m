
paramSet = {[0.03   0.00085  250]};
%  smallest grid size epsilon  max wavenumber (just keep it at like 150)
 
   parami = paramSet{1};
   
fd=@(p) drectangle(p,-2, 2, -0.5,0.5);
fh=@(p) min(parami(1)+max(0,0.5*drectangle(p, -2, 0, -0.5, 0).^3),0.1);
 %fh=@(p) min(max(parami(1),parami(1)+1.4*drectangle(p, 1.99,2.01, 0.05, 0.5).^3),0.1);
%fh=@(p) parami(1);
[p,t]=distmesh2d(fd,fh,parami(1),[-2, -0.5; 2, 0.5],[-2  -0.5; -2  0.5; 2  0.5;  2  -0.5]);
% this for cross section plot
% [X,Y] =    meshgrid(-0.39:0.02:-0.23);
% X = reshape(X,numel(X),1);
% Y = reshape(Y,numel(Y),1);
% xp = [X,Y]';
% xp = [[0.2,0.2];[0.4,0.3]];
[X,Y] =    meshgrid(-1.6:0.11:0,-0.5:0.06:0);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
[X2,Y2] =    meshgrid(-1.9:0.05:-1.6,-0.5:0.05:0);
X2 = reshape(X2,numel(X2),1);
Y2 = reshape(Y2,numel(Y2),1);
xp2 = [X2,Y2]';
xp=[xp,xp2];
fdtest=@(p) drectangle(p,-2, 0, -0.5, 0);
xp=xp(:,fdtest(xp')<-0.02); %mess here

velocities = zeros(length(paramSet), length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);

for i = 1:length(paramSet)
    parami = paramSet{i};
  

maxWaveNum = parami(3);
L = 12;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told]=velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L,1,2);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
%
[p,t]=distmesh2d(fd,@huniform,parami(1),[-2, -0.5; 2, 0.5],[-2  -0.5; -2  0.5; 2  0.5;  2  -0.5]);
figure
hold on
v = [0  0 ; 0  1; 1   1; 1  0];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none');

quiver(xp(1,:),xp(2,:),velocities(1,:,1),velocities(1,:,2),'b')
%quiver(xp(1,:),xp(1,:),-velocities(1,:,2),velocities(1,:,1),'b')
fullXP = [xp(1, :)',xp(2, :)';
-xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)';
-xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
-fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2);
-fuck(:,1),-fuck(:,2)];
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.04;
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
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
finalIC = zeros(length(p),1);
finalIC(ix) = ic;
figure
plotFESol(p,t,finalIC)
hold on
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-2  -0.5 ; 2  -0.5; 2   0.5; 2  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)