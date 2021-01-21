function  triangles = nearbyTriangles(t,nodes)
triangles = [];
for i = 1:numel(nodes)
    triangles = [triangles;find(t(:,1) == nodes(i) | t(:,2) == nodes(i)...
                               | t(:,3) == nodes(i))];
end
triangles = unique(triangles);
end
