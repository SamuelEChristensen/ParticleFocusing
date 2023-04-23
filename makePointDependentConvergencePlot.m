figure
hold on
for i=2:5
    for j=1:4
        if i==5 && j==4
            continue
        end
        scatter(A{i,j}(:,1),A{i,j}(:,2),'MarkerFaceColor', [1  j/5  0],'MarkerEdgeColor', [1  j/5  0])
    end
end
axis image
title('focusing wrt sampling')

for i=1:5
    for j=1:4
        if i==5 && j==4
            VconvClose(i,j,1) = NaN;
            VconvClose(i,j,2) = NaN;
            VconvFar(i,j,1) = NaN;
            VconvFar(i,j,2) = NaN;
            continue
        end
        v = C{i,j};
        xp = B{i,j};
        Fy = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,2),'natural');
        Fx = scatteredInterpolant(xp(1,:)', xp(2,:)', v(:,1),'natural');
        fuckShit = @(t,y) [Fx(y(1,:),y(2,:));Fy(y(1,:),y(2,:))];
        VconvClose(i,j,:) = fuckShit(1,[0.2,-0.354]');
        VconvFar(i,j,:) = fuckShit(1,[0.15,-0.254]');
    end
end
figure
hold on
for i = 1:4
    plot(1:5, VconvClose(:,i,1))
end
title('convergence on slow manifold')
legend('sample spacing 0.1','sample spacing 0.08','sample spacing 0.06','sample spacing 0.04')
xlabel('increasing mesh density')
figure
hold on
for i = 1:4
    plot(1:5, VconvFar(:,i,1))
end
xlabel('increasing mesh density')
title('convergence on fast manifold')
legend('sample spacing 0.1','sample spacing 0.08','sample spacing 0.06','sample spacing 0.04','location','se')
