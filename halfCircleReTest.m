c=1.2;
ij=0;
for Re = [ 15 ]
    ij=ij+1;
    fd=@(p) ddiff(sqrt((p(:,1)/c).^2+p(:,2).^2)-1,drectangle(p,-1,0,-1,1));
    bbox = [0,-1;c,1];
    
    fh=@(p) min(0.035+max(0,0.5*drectangle(p, 0, c, 0.15, 0.5).^3), min(0.02+max(0,0.5*drectangle(p, 0, c, 0, 0.15).^3), 0.075));
    pfix = [0,-1;0,1];
    [X,Y] =    meshgrid(0:0.11:c,0.2:0.11:1);
    X = reshape(X,numel(X),1);
    Y = reshape(Y,numel(Y),1);
    xp = [X,Y]';
    [X2,Y2] =    meshgrid(0.03:0.12:c,0.03:0.8:0.2);
    X2 = reshape(X2,numel(X2),1);
    Y2 = reshape(Y2,numel(Y2),1);
    xp2 = [X2,Y2]';
    xp = [xp,xp2];
    xp = unique(xp','rows');
    xp=xp';
    
    settings = {'reflecty'};
    figure
    [xpi,vi, eqPointsi] = findMigrationCurveChannel(fd,bbox,pfix, Re,fh,xp,settings);
    EQPstore{ij} = eqPointsi;
    vstore{ij} = vi;
    xpStore{ij} = xpi;
end