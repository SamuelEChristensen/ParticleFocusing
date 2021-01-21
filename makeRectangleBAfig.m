%%%load square
load squareChannel.mat
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
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)

%%% load 1.5x1
load 1.5by1.mat
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
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-0.75  -0.5 ; 0.75  -0.5; 0.75   0.5; -0.75  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)




%%% load 2x1
load 2by1Channel.mat
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
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-1  -0.5 ; 1  -0.5; 1   0.5; -1  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)






%%% load 2.5x1
load 2.5by1Channel.mat
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
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',1.5);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)

% %%% load 3x1
% load 3by1by10Channel.mat
% figure
% hold on
% v = [0  0 ; 0  1; 1   1; 1  0];
% patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none');
% 
% quiver(xp(1,:),xp(2,:),velocities(1,:,1),velocities(1,:,2),'b')
% %quiver(xp(1,:),xp(1,:),-velocities(1,:,2),velocities(1,:,1),'b')
% fullXP = [xp(1, :)',xp(2, :)';
% -xp(1, :)',xp(2, :)';
% xp(1, :)',-xp(2, :)';
% -xp(1, :)',-xp(2, :)'];
% fuck = reshape(velocities(1,:,:),length(xp),2);
% fullVelocities = [fuck(:,1),fuck(:,2);
% -fuck(:,1),fuck(:,2);
% fuck(:,1),-fuck(:,2);
% -fuck(:,1),-fuck(:,2)];
% Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
% Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
% fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
% ix = fd(p)<-0.04;
% sp=p(ix,:);
% output = zeros(size(sp));
% finalAns = zeros(size(sp));
% figure
% hold on
% for i = 1:length(sp)
% [T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,300], sp(i,:));
% finalAns(i,:) = Y(end,:);
% plot(Y(:,1),Y(:,2))
% end
% [eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
% finalIC = zeros(length(p),1);
% finalIC(ix) = ic;
% figure
% plotFESol(p,t,finalIC)
% hold on
% scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
% v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
% patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
% axis image
% set(gca,'xtick',[],'ytick',[])
% xlabel('x','fontSize',24)
% ylabel('y','fontSize',24)



%%% load 4x1
load 4by1Channel.mat
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
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)