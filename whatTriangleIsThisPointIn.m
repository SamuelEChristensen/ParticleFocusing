function n = whatTriangleIsThisPointIn(p, t, x)
n = 0;
    for i = 1:length(t)
        nodes = t(i, 1:3);
        P = [ones(1, 3); p(:, nodes)];
        areaOfElement = abs(det(P)) / 2;
        
        P1 = [ones(1, 3); p(:, nodes([1,2]))  x];
        areaOfP1 = abs(det(P1)) / 2;
        
        P2 = [ones(1, 3); p(:, nodes([1,3]))  x];
        areaOfP2 = abs(det(P2)) / 2;
        
        P3 = [ones(1, 3); p(:, nodes([3,2]))  x];
        areaOfP3 = abs(det(P3)) / 2;
        if abs(areaOfP1+areaOfP2+areaOfP3-areaOfElement) < 10000*eps
            n=i;
            break
        end
            
    end
    if n==0
        fprintf('The point was not in an element')
    end
end