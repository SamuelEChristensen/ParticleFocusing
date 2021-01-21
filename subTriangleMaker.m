function [subp,subt] = subTriangleMaker(points)

mids = [mean(points([1,2],:),1)
        mean(points([1,3],:),1)
        mean(points([3,2],:),1)];


subp = [points
        mids
        ];
subt = [1  4  5
        2  4  6
        3  5  6
        4  5  6];
end