clear EQPstore vstore xpStore
pv = [-0.5  -0.5; -0.5  0.5; 0.5  0.5;  0.5  -0.5;  -0.5  -0.5];
fh=@(p) min(0.025+max(0,0.5*drectangle(p, -2, 0, -0.5, 0).^3),0.05);  % custom mesh distance that puts more points in one quadrant

settings = {'reflectx','reflecty'};  %commands to take advantage of symmetry
Re = 1;
count1=0;
count2=0;
for hi = [0.1  0.08  0.06  0.04]
    for pi = [0.03]
        if hi==0.04 && pi == 0.015
            continue
        end
        [X,Y] =    meshgrid(-0.5:hi:-0.01,-0.5:hi:-0.01);
        X = reshape(X,numel(X),1);
        Y = reshape(Y,numel(Y),1);
        xp = [X,Y]';  %particle testing locations
        fh=@(p) min(pi+max(0,0.5*drectangle(p, -2, 0, -0.5, 0).^3),0.075);  % custom mesh distance that puts more points in one quadrant
        count1=count1+1;
        count2=count2+1;
        figure
        [xpi,vi, eqPointsi] = findMigrationPolyChannel(pv, Re, fh,xp, settings);
        EQPstore{count1,count2} = eqPointsi;
        vstore{count1,count2} = vi;
        xpStore{count1,count2} = xpi;
    end
end
figure
hold on
for i=1:4
    scatter(EQPstore{i,i}(:,1),EQPstore{i,i}(:,2))
end
% 
% ix = fd(p)<-0.1;
% pix = p(ix,:);
% VVar = zeros(length(pix),1);
% fullVVar = zeros(length(p),1);
% Vguess = zeros(5,4,2);
% for h=1:length(pix)
%     for i=1:5
%         for j=1:4
%             if i<4 || j<3
%                 Vguess(i,j,1) = NaN;
%                 Vguess(i,j,2) = NaN;
%                 continue
%             end
%             if i==5 && j==4
%                 Vguess(i,j,1) = NaN;
%                 Vguess(i,j,2) = NaN;
%                 continue
%             end
%             v = C{i,j};
%             xp = B{i,j};
%             Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
%             Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
%             fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
%             Vguess(i,j,:) = fuckShit(1,pix(h,:)');
%         end
%     end
%     VVar(h) = min((var(Vguess(:,:,1),0,'all', 'omitnan') + var(Vguess(:,:,2),0, 'all','omitnan')).^0.5/(mean((Vguess(:,:,1).^2+Vguess(:,:,2).^2).^0.5,'all','omitnan')),2);
% end
% fullVVar(ix) = VVar;
% figure
% plotFESol(p,t,fullVVar)