L = 5;
xp = [0; 0];  %particle testing locations
Re = 1;
gamma = -0*3.8;
delta = -8;
count=0;
epsilon=0.00085;
maxWaveNum=150;
for hi=1:10
    pv = 4*[-0.5  -0.5; -0.5  0.5; 0.5-0.05*(hi-1)  0.5;  0.5-0.05*(hi-1)  -0.5;  -0.5  -0.5];
    bbox=[min(pv(:,1)),min(pv(:,2));
        max(pv(:,1)), max(pv(:,2))];
    fd=@(p) dpoly(p,pv);
    xlength = bbox(2)-bbox(1);
    ylength = bbox(4)-bbox(3);
    fh=@(p) min(min(0.01+max(0,1/3*drectangle(p, -0.001, 0.01, -0.001, 0.001).^2),0.08 + max(0,1/20*drectangle(p, -0.5, 0.5, -0.5, 0.5))),0.3);  % custom mesh distance that puts more points in one quadrant
    smallScale = min(fh(xp'));
    disp(['smallScale = ',num2str(smallScale)])
    sol=@(x) zeros(size(x(:,1)));
    f=@(x) -ones(size(x(:,1)));
    fb=@(p) sol(p);
    
    [p,t]=distmesh2d(fd,fh, smallScale,bbox,pv);
    count=count+1;
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    [velocity,~,~] = velocitySolveRaymond(p, t, epsilon, waveNumbers, xp, L, gamma, delta, Re);
    vstore{count} = (2*pi)^0.5*maxWaveNum/L*2*velocity*Re;
    ray1 = 0.2177*abs(gamma)*abs(delta)^(2/3)*Re^(2/3);
    ray2 = ray1 + 0.0031*gamma^2*Re;
    ray3 = ray2 - 0.0322*gamma^3*Re^(4/3)*abs(delta)^(-2/3);
end
