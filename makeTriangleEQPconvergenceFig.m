pv = [-0.5  0;  0  (3/4)^0.5;  0.5  0; -0.5  0];
Re = 1;
fh=@(p) min(0.025+max(0,0.5*drectangle(p, -2, -0.1, -0, 1).^3),0.125);  % custom mesh distance that puts more points in one quadrant

xpStore=cell(6,1);
vStore = cell(6,1);
eqpStore=cell(6,1);
count=1;
for hi=[0.08,0.07,0.06,0.05,0.04,0.032, 0.025] %0.08,0.07,0.06,0.05,0.04,
[X,Y] =    meshgrid(-0.75:hi:0,0.02:hi:1);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);
xp = [X,Y]';  %particle testing locations
settings = {'reflectx'};  %commands to take advantage of symmetry
[xpi,vi, eqpi] = findMigrationPolyChannel(pv, Re, fh, xp, settings);
xpStore{count} = xpi;
vStore{count} = vi;
eqpStore{count} = eqpi;
count=count+1;
end

error = zeros(1,8);
figure
hold on
for i=1:7
scatter(eqpStore{i}(:,1), eqpStore{i}(:,2))
error(i) = norm(eqpStore{i}-eqpStore{7})/(norm(eqpStore{7}));
end
figure
plot((log([0.09, 0.08,0.07,0.06,0.05,0.04,0.032])/log(2)), log2(error(1:7)),'lineWidth',2.0)
xlabel('Grid Spacing')
ylabel('Rel. Error of Focusing Position')
XTickLabels = flipud(cellstr(num2str(round(10*log2([0.09, 0.08,0.07,0.06,0.05,0.04,0.032])')/10, '2^{%d}')));
YTickLabels = (cellstr(num2str([-3.5,-3,-2.5,-2]', '2^{%d}')));
ax = gca;
ax.FontSize = 13;
set(gca,'xticklabel',XTickLabels,'xtickmode','manual');
set(gca,'yticklabel',YTickLabels,'ytickmode','manual');
xticks(fliplr([2.^round(log2([0.09, 0.08,0.07,0.06,0.05,0.04,0.032]))]))
yticks([2^(-3.5),2^(-3),2^(-2.5),2^(-2)])