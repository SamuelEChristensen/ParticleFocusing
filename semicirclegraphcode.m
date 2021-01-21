c=1.2;
 bbox = [0,-1;c,1];
pfix = [0,-1;0,1];
fd=@(p) ddiff(sqrt((p(:,1)/c).^2+p(:,2).^2)-1,drectangle(p,-2,0,-1,1));
[p,t]=distmesh2d(fd,@huniform,0.03*max(2,1.2) ,bbox,pfix);
xp = 0.96*[[c*cos(0.05:0.1:2*pi);sin(0.05:0.1:2*pi)],0.9*[c*cos(0.05:0.1:2*pi);sin(0.05:0.1:2*pi)],0.8*[c*cos(0.05:0.1:2*pi);sin(0.05:0.1:2*pi)],0.7*[c*cos(0.05:0.2:2*pi);0.4*sin(0.05:0.2:2*pi)],0.6*[c*cos(0.05:0.3:2*pi);sin(0.05:0.3:2*pi)],0.5*[c*cos(0.05:0.4:2*pi);sin(0.05:0.4:2*pi)],0.4*[c*cos(0:0.5:2*pi);sin(0:0.5:2*pi)],0.3*[c*cos(0:0.6:2*pi);sin(0:0.6:2*pi)],0.2*[c*cos(0:0.7:2*pi);sin(0:0.7:2*pi)],0.1*[c*cos(0:0.8:2*pi);sin(0:0.8:2*pi)]]
fullFigXP = xp(:,fd(xp')<-0.03);
fullXP=xpStore{4}';
v = vstore{4};
fullVelocities = v;
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.02;
sp=fullFigXP';
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,200], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2),'b','LineWidth',2)
arowscal = 10;
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
%finalIC = zeros(length(sp),1);
%finalIC(ix) = ic;
%figure
%plotFESol(p,t,finalIC)
hold on
scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),60,'b','filled')
%set(gca,'xtick',[],'ytick',[])
ax = gca;
ax.FontSize = 16; 


xlabel('x','fontSize',24)
ylabel('y','fontSize',24)



c=1;
fd=@(p) ddiff(sqrt((p(:,1)/c).^2+p(:,2).^2)-1,drectangle(p,-1,0,-1,1));
xp = 0.97*[[c*cos(0:0.1:2*pi);sin(0:0.1:2*pi)],0.9*[c*cos(0:0.1:2*pi);sin(0:0.1:2*pi)],0.8*[c*cos(0:0.1:2*pi);sin(0:0.1:2*pi)],0.7*[c*cos(0:0.2:2*pi);0.4*sin(0:0.2:2*pi)],0.6*[c*cos(0:0.3:2*pi);sin(0:0.3:2*pi)],0.5*[c*cos(0:0.4:2*pi);sin(0:0.4:2*pi)],0.4*[c*cos(0:0.5:2*pi);sin(0:0.5:2*pi)],0.3*[c*cos(0:0.6:2*pi);sin(0:0.6:2*pi)],0.2*[c*cos(0:0.7:2*pi);sin(0:0.7:2*pi)],0.1*[c*cos(0:0.8:2*pi);sin(0:0.8:2*pi)]]
fullFigXP = xp(:,fd(xp')<-0.02);
fullXP=xpStore{1}';
v = vstore{1};
fullVelocities = v;
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.04;
sp=fullFigXP';
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,200], sp(i,:));
finalAns(i,:) = Y(end,:);
plot(Y(:,1),Y(:,2),'r','LineWidth',2)
arowscal = 10;
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
%finalIC = zeros(length(sp),1);
%finalIC(ix) = ic;
%figure
%plotFESol(p,t,finalIC)
scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),60,'r','filled')
%set(gca,'xtick',[],'ytick',[])

ax = gca;
ax.FontSize = 16; 
xlabel('x','fontSize',24)
ylabel('y','fontSize',24)