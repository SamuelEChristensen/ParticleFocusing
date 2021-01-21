figXP = [[0.9+0*(0:0.14:0.8);0:0.14:0.8], [0:0.14:0.8;0.9+0*(0:0.14:0.8)], [0.1+0*(0:0.14:1);(0:0.14:1)], [(0:0.14:1);0.1+0*(0:0.14:1)],[0.14;0.15], [0.14;0.13]]/2;
fullFigXP = [figXP(1, :)',figXP(2, :)';
-figXP(1, :)',figXP(2, :)';
figXP(1, :)',-figXP(2, :)';
-figXP(1, :)',-figXP(2, :)'];
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
sp=fullFigXP;
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,100], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2),'b')
arowscal = 10000;
quiver(Y(:,1), Y(:,2), Fx(Y), Fy(Y),'b','LineWidth',1.5)
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
%finalIC = zeros(length(sp),1);
%finalIC(ix) = ic;
%figure
%plotFESol(p,t,finalIC)
hold on
scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),110,'k','filled')
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)