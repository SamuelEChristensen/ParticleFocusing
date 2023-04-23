%clear EQPstore vstore xpStore
pv = [-0.5  -0.5; -0.5  0.5; 0.5  0.5;  0.5  -0.5;  -0.5  -0.5];

xp = [-0.2; -0.1];  %particle testing locations
settings = {'reflectx','plotOff'};  %commands to take advantage of symmetry
Re = 1;
count=0     ;
for Re = [1,50]
    pi=[0.005];
    fh=@(p) min(pi+max(0,0.5*drectangle(p, -0.19, 0.-0.21, -0.09, -0.11).^3),0.08);  % custom mesh distance that puts more points in one quadrant

    count=count+1;
    [xpi,vi, eqPointsi] = findMigrationPolyChannel(pv, Re, fh, xp, settings);
    EQPstore{count} = eqPointsi;
    vstore{count} = vi;
    xpStore{count} = xpi;
end
figure
hold on
plot(log(1:512)/log(2),log(((vstore{1}(1:512,1)-vstore{1}(1024,1)).^2+(vstore{1}(1:512,2)-vstore{1}(1024,2)).^2).^0.5/norm(vstore{1}(1024,:)))/log(2),'lineWidth',2.0)
%plot(log(1:512),log(((vstore{2}(1:512,1)-vstore{2}(1024,1)).^2+(vstore{2}(1:512,2)-vstore{2}(1024,2)).^2).^0.5))
plot(log(1:512)/log(2),log(((vstore{2}(1:512,1)-vstore{2}(1024,1)).^2+(vstore{2}(1:512,2)-vstore{2}(1024,2)).^2).^0.5/norm(vstore{2}(1024,:)))/log(2),'lineWidth',2.0)
ax = gca;
ax.FontSize = 13;
XTickLabels = (cellstr(num2str([0,3,6,9]', '2^{%d}')));
YTickLabels = (cellstr(num2str([-10,-7,-4,-1]', '2^{%d}')));
set(gca,'xticklabel',XTickLabels,'xtickmode','manual');
set(gca,'yticklabel',YTickLabels,'ytickmode','manual');
set(gca, 'box','on')
xticks([[0,(3),(6),9]])
yticks([(-10),(-7),(-4),(-1)])
xlabel('Number Of Fourier Modes Used','fontSize',18)
ylabel('Relative Error','fontSize',18)