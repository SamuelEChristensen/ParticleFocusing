pastel = [154, 147, 255];
figure
hold on
xvec = reshape(x, numel(x),1);
yvec = reshape(y, numel(y),1);
fullX = [xvec;
-xvec;
xvec;
-xvec];

fullY = [yvec;
yvec;
-yvec;
-yvec];


fuckvec = [reshape(uvec,numel(x),1),reshape(vvec,numel(x),1)];

fullvec = [fuckvec(:,1),fuckvec(:,2);
-fuckvec(:,1),fuckvec(:,2);
fuckvec(:,1),-fuckvec(:,2);
-fuckvec(:,1),-fuckvec(:,2)];


fuckf = [uf', vf'];
fullf = [fuckf(:,1),fuckf(:,2);
-fuckf(:,1),fuckf(:,2);
fuckf(:,1),-fuckf(:,2);
-fuckf(:,1),-fuckf(:,2)];
figure 
hold on
%contour(fullX, fullY,(fuckvec(:,1).^2+fuckvec(:,2).^2).^0.5,'k')
%contour(x,y,(reshape(uf,9,9).^2+reshape(vf,9,9).^2).^0.5,'color',pastel/255)
pastel = [154, 147, 255];
quiver(fullX,fullY,fullf(:,1),fullf(:,2),'color','k','lineWidth',3.0)
quiver(fullX,fullY,fullvec(:,1),fullvec(:,2),'color', pastel/255,'lineWidth',3.0)
v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',2.0);
set(gca,'xtick',[],'ytick',[])
xlabel('x')
ylabel('y')
