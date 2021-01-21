clear EQPstore vstore xpStore
pv = [-0.5  -0.5; -0.5  0.5; 0.5  0.5;  0.5  -0.5;  -0.5  -0.5];
fh=@(p) min(0.025+max(0,0.5*drectangle(p, -2, 0, -0.5, 0).^3),0.05);  % custom mesh distance that puts more points in one quadrant
[X,Y] =    meshgrid(-0.5:0.05:0,-0.5:0.05:0);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';  %particle testing locations
settings = {'reflectx','reflecty'};  %commands to take advantage of symmetry

count=0;
for Re = [1   10  40  70  100]
    count=count+1;
    [xpi,vi, eqPointsi] = findMigrationPolyChannel(pv, Re, fh, xp, settings);
    EQPstore{count} = eqPointsi;
    vstore{count} = vi;
    xpStore{count} = xpi;
end
figure
hold on
for i =1:5
    scatter(EQPstore{i}(:,1),EQPstore{i}(:,2))
end