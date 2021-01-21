 %pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;  %polygon vertices
%1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
t = linspace(0, 2*pi, 6)';
pv =  [cos(t),sin(t)];  %polygon vertices
  
Re = 1; %Reynolds number
[xp,v] = findMigrationPolyChannel(pv,1);

figure
hold on
patch('xData', pv(:,1),'ydata', pv(:,2),'faceColor','none');
quiver(xp(1,:),xp(2,:),v(:,1)',v(:,2)','b')
[p,t]=distmesh2d(fd,@huniform,0.03 ,bbox,pv);
Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
fd=@(p) dpoly(p,pv);
ix = fd(p)<-0.08;
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
patch('xData', pv(:,1),'ydata',pv(:,2),'faceColor','none','lineWidth',2.0);
axis image
set(gca,'xtick',[],'ytick',[])
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)