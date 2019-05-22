function plotStokesSol(p,t,Uwn, varargin)
% if size(p,1)~=length(p)
% p=p';
% end
% if size(U,1)~=length(U)
% U=U';
% end
zwn=length(Uwn)/2;
if nargin < 3 
    error('not enough inputs')
    
end
if nargin == 3 
    wn = zwn;
end
if nargin == 4
    wn = zwn + varargin{1};
end
if nargin > 4 
    wn = zwn;
    disp('too many inputs, using k=0')
end

    titles={'U1','U2','U3','P'};
figure
for i=1:4
subplot(2,2,i)
set(gca, 'Projection', 'perspective')
caxis([0 1])
daspect([1 1 2])
view(70, 20)
axis vis3d
axis([-1 1 -1 1 -1 1])
cla
trisurf(t, p(:,1), p(:,2), real(Uwn{wn,i}), 'EdgeColor', 'none', 'FaceColor', 'interp');
colorbar
title(titles{i})
drawnow
xlabel('X')
ylabel('Y')
end
end