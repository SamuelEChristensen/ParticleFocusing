
paramSet = {[0.014   0.001  250]};
%  smallest grid size epsilon  max wavenumber (just keep it at like 300)
 

% this for cross section plot
%[X,Y] =  meshgrid(-0.9:0.05:5);
[X,Y] =  meshgrid(-0.9:0.03:5,-0.9:0.06:5);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
xp=xp(:,fd(xp')<-0.01); %mess here
[X,Y] =  meshgrid(0.65:0.02:0.9,0.2:0.035:0.9);
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
L = 6;

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,pold,told]=velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',4,1,2);
velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
end
% 
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = fuck;
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
[X,Y] = meshgrid(0:0.04:1,0:0.04:1);
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(X),1);
quivPoints = [X,Y];
ix=fd(quivPoints)<-0.04;
figure
quiver(quivPoints(ix,1),quivPoints(ix,2),Fx(quivPoints(ix,:)),Fy(quivPoints(ix,:)),'k','lineWidth',1.0)
ip = fd(p)<-0.04;
im = fd(p)>-0.06;
iFinal = ip&im;
ip = fd(p)<-0.15;
im = fd(p)>-0.18;
iFinal = iFinal|(ip&im);
figXP = uniquetol(p(iFinal,:),0.04,'byRows',1);
sp=[figXP;[0.8,0.5]];
output = zeros(size(sp));
finalAns = zeros(size(sp));
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,100], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2),'b','lineWidth',2.0)
arowscal = 10000;
%quiver(Y(:,1), Y(:,2), Fx(Y)./(Fx(Y).^2+Fy(Y).^2), Fy(Y)./(Fx(Y).^2+Fy(Y).^2),'b','LineWidth',1.5)
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
%finalIC = zeros(length(sp),1);
%finalIC(ix) = ic;
%figure
%plotFESol(p,t,finalIC)
scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),50,'k','filled')
axis image
set(gca,'xtick',[],'ytick',[])
plot(wah(:,1)/1787,wah(:,2)/1787,'k','lineWidth',1.5)