function [xp, v, eqPoints] = findMigrationPolyChannel(pv, Re, fhUser, xpUser, settings, varargin)
%% Calculation
narginchk(2,5)
bbox=[min(pv(:,1)),min(pv(:,2));
    max(pv(:,1)), max(pv(:,2))];


paramSet = [0.02   0.00085  128];
%  smallest grid size epsilon  max wavenumber

tic
fd=@(p) dpoly(p,pv);
xlength = bbox(2)-bbox(1);
ylength = bbox(4)-bbox(3);
if nargin==2
    
    [p,t]=distmesh2d(fd,@huniform,0.0251*min(xlength,ylength),bbox,pv);
              
    [X,Y] =    meshgrid(bbox(1):0.05*xlength:bbox(2),bbox(3):0.05*ylength:bbox(4));
    X = reshape(X,numel(X),1);
    Y = reshape(Y,numel(Y),1);
    xp = [X,Y]';
    xp=xp(:,fd(xp')<-0.02); %Mess around here for boundary distance
    settings ={ };
end

if nargin>2
    smallScale = min(fhUser(xpUser'));
    disp(['smallScale = ',num2str(smallScale)])
    [p,t]=distmesh2d(fd,fhUser, smallScale,bbox,pv);
    xp=xpUser(:,fd(xpUser')<-0.02); %Mess around here for boundary distance
    if nargin<5
        settings ={ };
    end
end
toc
velocities = zeros(size(xp,2),2);
%velocities = zeros(paramSet(3),2);
sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));
fb=@(p) sol(p);




maxWaveNum = paramSet(3);
%L = 5*max(xlength,ylength);
L = 3*max(xlength,ylength);
waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
tic
[bgFlow,~,~] = poissonSolver(p,t,f,fb);
[u,~,~]=velocitySolveDiscDecomp(p,t, paramSet(2), waveNumbers, xp, bgFlow',L,1,2, Re);
%[u,~,~]=DisturbanceFieldSolveDiscDecomp(p,t, paramSet(2), waveNumbers, xp, bgFlow',L,1,2, Re);
velocities(:,:) = (2*pi)^0.5*maxWaveNum/L*u*2;
toc
v = velocities;
%% Plotting
if ismember('plotOff',settings)
    eqPoints=[];
end
if ~ismember('plotOff',settings)
    
    
    
    if ismember('reflectx',settings)
        v = [v; -v(:,1), v(:,2)];
        xp = [xp, [-xp(1,:); xp(2,:)]];
    end
    
    if ismember('reflecty',settings)
        v = [v; v(:,1), -v(:,2)];
        xp = [xp, [xp(1,:); -xp(2,:)]];
    end
    
    
    figure
    hold on
    patch('xData', pv(:,1),'ydata', pv(:,2),'faceColor','none');
    quiver(xp(1,:),xp(2,:),v(:,1)',v(:,2)','b')
    [p,t]=distmesh2d(fd,@huniform,0.029 ,bbox,pv);
    Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
    Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
    fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
    fd=@(p) dpoly(p,pv);
    ix = fd(p)<-0.03*max(xlength,ylength);
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
    plotFESol(p,t,finalIC)
    hold on
    scatter3(finalAns(:,1),finalAns(:,2),repmat(max(ic)+1,length(finalAns),1),110,'k','filled')
    patch('xData', pv(:,1),'ydata',pv(:,2),'faceColor','none','lineWidth',2.0);
    axis image
    set(gca,'xtick',[],'ytick',[])
end
end