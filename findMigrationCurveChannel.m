function [xp, v, eqPoints] = findMigrationCurveChannel(fd, bbox, pfix, Re, fhUser, xpUser, settings, varargin)
%% Calculation
narginchk(3,7)

paramSet = [0.04   0.00085  150];
%  smallest grid size epsilon  max wavenumber

xlength = bbox(2)-bbox(1);
ylength = bbox(4)-bbox(3);
if nargin==2
    
    [p,t]=distmesh2d(fd,@huniform,0.025*min(xlength,ylength),bbox);
    
    [X,Y] =    meshgrid(bbox(1):0.04*xlength:bbox(2),bbox(3):0.04*ylength:bbox(4));
    X = reshape(X,numel(X),1);
    Y = reshape(Y,numel(Y),1);
    xp = [X,Y]';
    xp=xp(:,fd(xp')<-0.02); %Mess around here for boundary distance
    settings ={ };
end

if nargin>2
    [p,t]=distmesh2d(fd,fhUser,0.02 ,bbox, pfix);
    xp=xpUser(:,fd(xpUser')<-0.02); %Mess around here for boundary distance
    if nargin<5
        settings ={ };
    end
end

velocities = zeros(length(xp),2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));
fb=@(p) sol(p);




maxWaveNum = paramSet(3);
L = 3*max(xlength,ylength);

waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,~,~]=velocitySolveDiscDecomp(p,t, paramSet(2), waveNumbers, xp, bgFlow',L,1,2, Re);
velocities(:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
%% Plotting
v = velocities;

if ismember('reflectx',settings)
    v = [v; -v(:,1), v(:,2)];
    xp = [xp, [-xp(1,:); xp(2,:)]];
end

if ismember('reflecty',settings)
    v = [v; v(:,1), -v(:,2)];
    xp = [xp, [xp(1,:); -xp(2,:)]];
end



[p,t]=distmesh2d(fd,@huniform,0.03*max(xlength,ylength) ,bbox,pfix);
Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
ix = fd(p)<-0.03*max(xlength,ylength);
sp=p(ix,:);
output = zeros(size(sp));
finalAns = zeros(size(sp));
%figure
%hold on
for i = 1:length(sp)
    [T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,600], sp(i,:));
    finalAns(i,:) = Y(end,:);
    %plot(Y(:,1),Y(:,2))
end
[eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1,'DataScale',1);
finalIC = zeros(length(p),1);
finalIC(ix) = ic;
figure
plotFESol(p,t,finalIC)
hold on
scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
axis image
set(gca,'xtick',[],'ytick',[])
end