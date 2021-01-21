%%%load square
load squareChannel.mat
fullXP = [xp(1, :)',xp(2, :)';
-xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)';
-xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
-fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2);
-fuck(:,1),-fuck(:,2)];
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];

sp = [0.025;0.025];


for i = 1:300
[T11,Y11] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y11(end,1))>0.1
    sp(1) = sp(1) -0.01>1/i;
end
 if abs(Y11(end,1))<0.1
     sp(1) = sp(1)  + 0.01>1/i;
 end   

 
end
figure
hold on
velocNorms = sum(fuckShit(0,Y11(T11<10,:)').^2,1).^0.5;
lastIndex11 = find(velocNorms==min(velocNorms));
plot(Y11(1:lastIndex11,1),Y11(1:lastIndex11,2),'b')

sp = [0.9,;0.9]/2;
for i = 1:300
[T12,Y12] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y12(end,1))>0.1
    sp(1) = sp(1) -0.025>1/i;
end
 if abs(Y12(end,1))<0.1
     sp(1) = sp(1)  + 0.025>1/i;
 end   

 
end
velocNorms = sum(fuckShit(0,Y12(T12<15,:)').^2,1).^0.5;
lastIndex12 = find(velocNorms==min(velocNorms));
plot(Y12(1:lastIndex12,1),Y12(1:lastIndex12,2),'b')


%%% load 1.5x1
load 1.5by1.mat
fullXP = [xp(1, :)',xp(2, :)';
-xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)';
-xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
-fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2);
-fuck(:,1),-fuck(:,2)];
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
sp = [0.2;0.025];
for i = 1:500
[T31,Y31] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y31(end,1))>0.1
sp(2) = sp(2) +0.01>1/i;
end
if abs(Y31(end,1))<0.1
sp(2) = sp(2)  - 0.01>1/i;                                      
end
end
figure
hold on
velocNorms = sum(fuckShit(0,Y31(T31<10,:)').^2,1).^0.5;
lastIndex31 = find(velocNorms==min(velocNorms));
plot(Y31(1:lastIndex31,1),Y31(1:lastIndex31,2),'b')
sp = [1.4,;0.9]/2;
for i = 1:500
[T32,Y32] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y32(end,1))>0.1
sp(2) = sp(2) +0.025>1/i;
end
if abs(Y32(end,1))<0.1
sp(2) = sp(2)  - 0.025>1/i;
end
end
velocNorms = sum(fuckShit(0,Y32(T32<10,:)').^2,1).^0.5;
lastIndex32 = find(velocNorms==min(velocNorms));
plot(Y32(:,1),Y32(:,2),'b')








%%% load 2x1
load 2by1Channel.mat
fullXP = [xp(1, :)',xp(2, :)';
-xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)';
-xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
-fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2);
-fuck(:,1),-fuck(:,2)];
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];

sp = [0.425;0.025];


for i = 1:300
[T21,Y21] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y21(end,1))>0.1
    sp(2) = sp(2) +0.01>1/i;
end
 if abs(Y21(end,1))<0.1
     sp(2) = sp(2)  - 0.01>1/i;
 end   

 
end
figure
hold on
velocNorms = sum(fuckShit(0,Y21(T21<20,:)').^2,1).^0.5;
lastIndex21 = find(velocNorms==min(velocNorms));
plot(Y21(1:lastIndex21,1),Y21(1:lastIndex21,2),'b')

sp = [1.9,;0.9]/2;
for i = 1:300
[T22,Y22] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y22(end,1))>0.1
    sp(2) = sp(2) +0.025>1/i;
end
 if abs(Y22(end,1))<0.1
     sp(2) = sp(2)  - 0.025>1/i;
 end   

 
end
velocNorms = sum(fuckShit(0,Y22(T22<20,:)').^2,1).^0.5;
lastIndex22 = find(velocNorms==min(velocNorms));
plot(Y22(1:lastIndex22,1),Y22(1:lastIndex22,2),'b')


% % %% load 2.5x1
% %  %% load 2.5x1
% % load 2.5by1Channel.mat
% % fullXP = [xp(1, :)',xp(2, :)';
% % -xp(1, :)',xp(2, :)';
% % xp(1, :)',-xp(2, :)';
% % -xp(1, :)',-xp(2, :)'];
% % fuck = reshape(velocities(1,:,:),length(xp),2);
% % fullVelocities = [fuck(:,1),fuck(:,2);
% % -fuck(:,1),fuck(:,2);
% % fuck(:,1),-fuck(:,2);
% % -fuck(:,1),-fuck(:,2)];
% % Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
% % Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
% % fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
% % sp = [0.75;0.025];
% % for i = 1:500
% % [T31,Y31] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
% % if abs(Y31(end,1))>0.1
% % sp(2) = sp(2) +0.1>1/i;
% % end
% % if abs(Y31(end,1))<0.1
% % sp(2) = sp(2)  - 0.1>1/i;
% % end
% % end
% % figure
% % hold on
% % velocNorms = sum(fuckShit(0,Y31(T31<10,:)').^2,1).^0.5;
% % lastIndex31 = find(velocNorms==min(velocNorms));
% % plot(Y31(1:lastIndex31,1),Y31(1:lastIndex31,2),'b')
% % 
% % sp = [2.4,;0.95]/2;
% % for i = 1:500
% % [T,Y32] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
% % if abs(Y32(end,1))>0.1
% % sp(2) = sp(2) -0.025>1/i;
% % end
% % if abs(Y32(end,1))<0.1
% % sp(2) = sp(2)  +  0.025>1/i;
% % end
% % end
% % 
% % plot(Y32(:,1),Y32(:,2),'b')

%%% load 4x1
load 4by1Channel.mat
fullXP = [xp(1, :)',xp(2, :)';
-xp(1, :)',xp(2, :)';
xp(1, :)',-xp(2, :)';
-xp(1, :)',-xp(2, :)'];
fuck = reshape(velocities(1,:,:),length(xp),2);
fullVelocities = [fuck(:,1),fuck(:,2);
-fuck(:,1),fuck(:,2);
fuck(:,1),-fuck(:,2);
-fuck(:,1),-fuck(:,2)];
Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];

sp = [1.425;0.025];

for i = 1:300
[T41,Y41] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y41(end,1))>1.4
    sp(2) = sp(2) +0.01>1/i;
end
 if abs(Y41(end,1))<1.4
     sp(2) = sp(2)  - 0.01>1/i;
 end   

 
end
figure
hold on
velocNorms = sum(fuckShit(0,Y41(T41<10,:)').^2,1).^0.5;
lastIndex41 = find(velocNorms==min(velocNorms));
plot(Y41(1:lastIndex41,1),Y41(1:lastIndex41,2),'b')

sp = [3.9,;0.9]/2;
for i = 1:300
[T42,Y42] = ode23s(@(t,y) fuckShit(t,y), [0,400], sp);
if abs(Y42(end,1))>1.4
    sp(2) = sp(2) +0.025>1/i;
end
 if abs(Y42(end,1))<1.4
     sp(2) = sp(2)  - 0.025>1/i;
 end   

 
end
velocNorms = sum(fuckShit(0,Y42(T42<10,:)').^2,1).^0.5;
lastIndex42 = find(velocNorms==min(velocNorms));
plot(Y42(1:lastIndex42,1),Y42(1:lastIndex42,2),'b')

figure
hold on
plot(-Y11(1:lastIndex11,1)+0.5,Y11(1:lastIndex11,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y12(1:lastIndex12,1)+0.5,Y12(1:lastIndex12,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y21(1:lastIndex21,1)+0.5+0.5,Y21(1:lastIndex21,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y22(1:lastIndex22,1)+0.5+0.5,Y22(1:lastIndex22,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y31(1:lastIndex31,1)+0.25+0.5,Y31(1:lastIndex31,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y32(1:lastIndex32,1)+0.25+0.5,Y32(1:lastIndex32,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y41(1:lastIndex41,1)+1.5+0.5,Y41(1:lastIndex41,2),'k','lineWidth',1.0,'MarkerSize',10)
plot(-Y42(1:lastIndex42,1)+1.5+0.5,Y42(1:lastIndex42,2),'k','lineWidth',1.0,'MarkerSize',10)
plot(-Y11(1:lastIndex11,1)+0.5, Y11(1:lastIndex11,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y12(1:lastIndex12,1)+0.5,Y12(1:lastIndex12,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y21(1:lastIndex21,1)+0.5+0.5,Y21(1:lastIndex21,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y22(1:lastIndex22,1)+0.5+0.5,Y22(1:lastIndex22,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y31(1:lastIndex31,1)+0.25+0.5,Y31(1:lastIndex31,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y32(1:lastIndex32,1)+0.25+0.5,Y32(1:lastIndex32,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y41(1:lastIndex41,1)+1.5+0.5,Y41(1:lastIndex41,2),'k','lineWidth',1.0)
plot(-Y42(1:lastIndex42,1)+1.5+0.5,Y42(1:lastIndex42,2),'k','lineWidth',1.0)
plot(-Y11(1:lastIndex11,1)+0.5,-Y11(1:lastIndex11,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y12(1:lastIndex12,1)+0.5,-Y12(1:lastIndex12,2),'b','lineStyle','--','lineWidth',1.0)
plot(-Y21(1:lastIndex21,1)+0.5+0.5,-Y21(1:lastIndex21,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y22(1:lastIndex22,1)+0.5+0.5,-Y22(1:lastIndex22,2),'r--o','lineWidth',1.0,'MarkerSize',5)
plot(-Y31(1:lastIndex31,1)+0.25+0.5,-Y31(1:lastIndex31,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y32(1:lastIndex32,1)+0.25+0.5,-Y32(1:lastIndex32,2),'g-->','lineWidth',1.0,'MarkerSize',5)
plot(-Y41(1:lastIndex41,1)+1.5+0.5,-Y41(1:lastIndex41,2),'k','lineWidth',1.0)
plot(-Y42(1:lastIndex42,1)+1.5+0.5,-Y42(1:lastIndex42,2),'k','lineWidth',1.0)


v = [-0.5  -0.5 ; 0.5  -0.5; 0.5   0.5; -0.5  0.5];
%patch('xData', v(:,1),'ydata',v(:,2),'faceColor','none','lineWidth',1.0);
plot(-v(:,1)+0.5,v(:,2),'k','lineWidth',2.0)
axis image
%set(gca,'xtick',[],'ytick',[])
xlabel('Distance From Channel Wall','fontSize',16)
ylabel('Height','fontSize',16)
% figXP = [[0.9+0>(0:0.15:1);0:0.15:1], [0:0.15:1;0.9+0>(0% figXP = [[0.9+0>(0:0.15:1);0:0.15:1], [0:0.15:1;0.9+0>(0:0.15:1)], [0.1+0>(0:0.15:1);(0:0.15:1)], [(0:0.15:1);0.1+0>(0:0.15:1)],[0.15;0.16], [0.15;0.14]]/2
% fullXP = [xp(1, :)',xp(2, :)';
% -xp(1, :)',xp(2, :)';
% xp(1, :)',-xp(2, :)';
% -xp(1, :)',-xp(2, :)'];
% fuck = reshape(velocities(1,:,:),length(xp),2);
% fullVelocities = [fuck(:,1),fuck(:,2);
% -fuck(:,1),fuck(:,2);
% fuck(:,1),-fuck(:,2);
% -fuck(:,1),-fuck(:,2)];
% Fy = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,2),'natural');
% Fx = scatteredInterpolant(fullXP(:,1),fullXP(:,2), fullVelocities(:,1),'natural');
% fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
% ix = fd(p)<-0.04;
% sp=figXP';
% output = zeros(size(sp));
% finalAns = zeros(size(sp));
% figure
% hold on
% for i = 1:length(sp)
% [T,Y] = ode23s(@(t,y) fuckShit(t,y), [0,100], sp(i,:));
% finalAns(i,:) = Y(end,:);
% plot(Y(:,1),Y(:,2),'b')
% arowscal = 10;
% quiver(Y(:,1), Y(:,2), Fx(Y), Fy(Y),'b')
% end
% [eqPoints,ia,ic] = uniquetol(finalAns, 0.1,'ByRows',1);
% %finalIC = zeros(length(sp),1);
% %finalIC(ix) = ic;
% %figure
% %plotFESol(p,t,finalIC)
% hold on
% scatter3(eqPoints(:,1),eqPoints(:,2),repmat(max(ic)+1,length(eqPoints),1),140,'k','filled')