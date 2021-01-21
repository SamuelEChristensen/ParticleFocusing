%test case for oseenSolver with gaussian approximation of delta function
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = {[0.01   0.01  512]  [0.0141   0.01  256]   [0.028   0.01  512]    [0.0141   0.01  512]};
paramSet = {[0.01   0.000085 150]};
%paramSet = {[0.005   0.002  1500] };
k=1:3;
pv = [cos(k*2*pi/3);sin(k*2*pi/3)]';
pv = [pv;pv(1,:)];
 fd =@(p) dpoly(p,pv);
    fh=@(p) min(parami(1)+max(0,0.3*drectangle(p, -1, 1, 0, 1).^3),0.05);
[X,Y] =  meshgrid(-1:0.08:1,0.02:0.08:1);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';
xp=xp(:,fd(xp')<-0.05);
velocities = zeros(length(paramSet), length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);
for i = 1:length(paramSet)
    parami = paramSet{i};
    
   
    fh=@(p) min(parami(1)+max(0,0.3*drectangle(p, -1, 1, 0, 1).^3),0.04);
    [p,t]= distmesh2d(fd,fh,0.04,[-1,-1; 1,1],pv);
    maxWaveNum = parami(3);
    L = 4;
    % z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
    % z = circshift(z, -maxWaveNum/2+1);
    %waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    tic
    [bgFlow,~,~] = poissonSolver(p,t,f,fb);
    [u,pold,told]=velocitySolveDiscDecomp(p,t, parami(2), waveNumbers, xp, bgFlow',L,1,2);
    velocities(i,:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
    toc
end
%
figure
x = -1:0.01:1;
hold on
plot(x,zeros(size(x)),'k')
plot(x,(1-x.^2).^0.5,'k')
quiver(xp(2,:),xp(1,:),velocities(1,:,2),velocities(1,:,1),'b')
quiver(-xp(2,:),xp(1,:),-velocities(1,:,2),velocities(1,:,1),'b')

[p,t]=distmesh2d(fd,@huniform,0.04,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
figure
hold on
v = [-1  -1 ; -1  1; 1   1; 1  -1];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none');

quiver(xp(1,:),xp(2,:),velocities(1,:,1),velocities(1,:,2),'b')
%quiver(xp(1,:),xp(1,:),-velocities(1,:,2),velocities(1,:,1),'b')
fullXP = [xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2)];
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
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1,'DataScale',1);
finalIC = zeros(length(p),1);
finalIC(ix) = ic;
figure
plotFESol([p(:,2),p(:,1)],t,finalIC)
hold on
scatter3(finalAns(:,2),finalAns(:,1),repmat(max(ic)+1,length(finalAns),1),140,'k','filled')
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)
axis equal