function plotFESol(p,t,U)
if size(p,1)~=length(p)
p=p';
end
if size(U,1)~=length(U)
U=U';
end
set(gca, 'Projection', 'perspective')
caxis([-1 1])
daspect([1 1 2])
axis vis3d
axis([-1 1 -1 1 -1 1])
cla
trisurf(t, p(:,1), p(:,2), real(U(1:length(p))), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'interp','LineStyle','None');
set(gca, 'Projection', 'Orthographic')
%colorbar
view(0, 90)
drawnow

end