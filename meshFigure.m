clear
fd=@(p) drectangle(p,-1,1,-1,1);
[p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
part = 1;
poop  = 1;
xpi = [0.0  -0.0]';
p=p';


numberOfNodes.old = size(p, 2);
numberOfElements = size(t, 1);
S = sparse(numberOfNodes.old, numberOfNodes.old);
counter = numberOfNodes.old + 1;
figure
for subfigIndex = 1:9
    poop  = floor((subfigIndex-1)/3);
    boop  = mod(subfigIndex-1,3)-1;
    for e = 1:numberOfElements
        nodes = t(e, :);
        if (S(nodes(1), nodes(2)) == 0)
            S(nodes(1), nodes(2)) = counter;
            S(nodes(2), nodes(1)) = counter;
            p(:, counter) = mean(p(:, [nodes(1) nodes(2)]), 2);
            counter = counter + 1;
        end
        if (S(nodes(2), nodes(3)) == 0)
            S(nodes(2), nodes(3)) = counter;
            S(nodes(3), nodes(2)) = counter;
            p(:, counter) = mean(p(:, [nodes(2) nodes(3)]), 2);
            counter = counter + 1;
        end
        if (S(nodes(1), nodes(3)) == 0)
            S(nodes(1), nodes(3)) = counter;
            S(nodes(3), nodes(1)) = counter;
            p(:, counter) = mean(p(:, [nodes(1) nodes(3)]), 2);
            counter = counter + 1;
        end
        t(e, 4) = S(nodes(1), nodes(2));
        t(e, 5) = S(nodes(2), nodes(3));
        t(e, 6) = S(nodes(1), nodes(3));
    end
    numberOfNodes.new = size(p, 2);
    
    
    xpn = whatTriangleIsThisPointIn(p,t,xpi);
    nearbyTris{part}=[xpn];
    for i =1:poop
        nearbyTris{part} = nearbyTriangles(t,t(nearbyTris{part},:));  %can be iterated for more surrounding layers
        if i== (poop-1)
            intNearbyTris{part} = nearbyTris{part};
        end
    end
    nodes = t(xpn, :);
    particleNodes{part} = nodes;
    
    
    if boop>-1
        for e = 1:numel(nearbyTris{part})
            [subp,subt] = subTriangleMaker(p(:,t(nearbyTris{part}(e),1:3))');
            allp=subp;
            allt=subt;
            for j = 1:boop
                allpi=[];
                allti=[];
                for i = 1:length(allt)  %iterate this part more for more subtriangles
                    [ssp,sst] = subTriangleMaker(allp(allt(i,:),:));
                    allpi =[allpi; ssp];
                    allti=[allti; sst+6*(i-1)];
                end
                allp = [allpi];
                allt = [allti];
            end
            
            supNodes = t(nearbyTris{part}(e), :);
            
            % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
            supP = [ones(1, 6);
                p(:, supNodes);
                p(1, supNodes).^2;
                p(1, supNodes) .* p(2, supNodes);
                p(2, supNodes).^2];
            
            allp=allp';
            F = size(allt,1);
            E = size(allp,2);
            [allp,allt] = makeMidPoints(allp,allt,F,E);
            tCell{e} = allt;
            pCell{e} = allp;
        end
    end
    
    figure
    X = p(1,:);
    Y = p(2, :);
    hold on
    scatter(X,Y,7*ones(size(X)),'k','filled')
    if boop>-1
        for i = 1:length(pCell)
            X = pCell{i}(1,:);
            Y = pCell{i}(2, :);
            scatter(X,Y,7*ones(size(X)),'k','filled')
        end
        for i = 1:length(tCell)
            for j = 1:size(tCell{i},1)
                
                X = pCell{i}(1,tCell{i}(j,1:3));
                Y = pCell{i}(2,tCell{i}(j,1:3));
                patch('xData',X,'yData', Y,'lineWidth',0.7,'faceColor','none')
            end
        end
    end
    scatter(X,Y,7*ones(size(X)),'k','filled')
    scatter(xpi(1),xpi(2),14,'r','filled')
    X = p(1,:);
    Y = p(2, :);
    for i = 1:length(t)
        patch('xData',X(t(i,1:3)),'yData', Y(t(i,1:3)),'lineWidth',0.7,'faceColor','none')
    end
    %title(['n=',num2str(boop+1),' D=',num2str(poop)])
    axis image
    set(gca,'XColor', 'none','YColor','none')
    nameStr = ['meshfigure',num2str(subfigIndex),'.png'];
    saveas(gcf, nameStr);
end