pv = [-0.5  -0.5; -0.5  0.5; 0.5  0.5;  0.5  -0.5;  -0.5  -0.5];
Re = 1;
fh=@(p) min(0.02+max(0,0.5*drectangle(p, -2, 0, -0.5, 0).^3),0.1);  % custom mesh distance that puts more points in one quadrant
[X,Y] =    meshgrid(-0.5:0.05:0,-0.5:0.05:0);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';  %particle testing locations
settings = {'reflectx','reflecty'};  %commands to take advantage of symmetry
[xp,v, eqPoints] = findMigrationPolyChannel(pv, Re, fh, xp, settings);
[X,Y] =    meshgrid(-0.5:0.05:0.5,-0.5:0.05:0.5);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp2 = [X,Y];  %particle testing locations
bbox=[min(pv(:,1)),min(pv(:,2));
max(pv(:,1)), max(pv(:,2))];
figure
hold on
patch('xData', pv(:,1),'ydata', pv(:,2),'faceColor','none');
quiver(xp(1,:),xp(2,:),v(:,1)',v(:,2)','b')
 xlength=1;
ylength=1;
Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
fd=@(p) dpoly(p,pv);
[p,t]=distmesh2d(fd,@huniform,0.055 ,bbox,pv);
ix = fd(xp2)<-0.07*max(xlength,ylength);
sp=xp2(ix,:);
output = zeros(size(sp));
finalAns = zeros(size(sp));
figure
hold on
for i = 1:length(sp)
[T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,600], sp(i,:));
finalAns(i,:) = Y(end,:);
z = zeros(size(Y(:,1)));
col = sum(fuckShit(0*Y',Y').^2,1).^0.5;  % This is the color, vary with x in this case.
surface([Y(:,1)';Y(:,1)'],[Y(:,2)';Y(:,2)'],[z';z'],[col;col],...
'facecol','no',...
'edgecol','interp',...
'linew',2);
end
patch('xData', pv(:,1),'ydata',pv(:,2),'faceColor','none','lineWidth',2.0);
axis image

%set(gca,'xtick',[],'ytick',[])
c = colorbar;
c.Label.String = 'Migration Speed';
c.Label.FontSize = 18;
scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(1)+1,length(eqPoints),1),80,'k','filled')