%reccomend doing on computer with lots of memory
% fd=@(p) drectangle(p,0,4,0,1);
% [p,t]=distmesh2d(fd,@huniform,0.02229,[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
% initialLength epsilon number of modes
%paramSet = { [0.0075   0.005  180]};
paramSet = {[0.009  0.0025  250]};
%xp = [2*ones(size(0.0:0.05:0.45)); 0.0:0.05:0.45];
xp = [2  ;0.25 ];
velocities = cell(length(paramSet), 1,2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);

for i = 1:length(paramSet)
    parami = paramSet{i};
    
    fd=@(p) drectangle(p,0,4,0,1);
    fh=@(p) min(max(parami(1), parami(1)+1.0*drectangle(p, 1.99, 2.01, 0.1, 0.5).^3),0.1);
    [p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
    maxWaveNum = parami(3);
    L = 4;
    % z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
    % z = circshift(z, -maxWaveNum/2+1);
    %waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    tic
    [bgFlow,~,~] = poissonSolver(p,t,f,fb);
    [u,p,t]=velocitySolveStressletDecomp(p,t, parami(2), waveNumbers, xp, bgFlow', 4);
    for j = 1:size(xp,2)
        velocities(i,j,1) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,1)*2};
        velocities(i,j,2) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,2)*2};
    end
    toc
    Ux = zeros(500,1);
    Uy = zeros(500,1);
    for j = 1:500
        xpi = [2; 1*j/501];
        xpn = whatTriangleIsThisPointIn(p,t,xpi);
        nodes = t(xpn, :);
        particleNodes{j} = nodes;
        % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
        P = [ones(1, 6);
            p(:, nodes);
            p(1, nodes).^2;
            p(1, nodes) .* p(2, nodes);
            p(2, nodes).^2];
        IPS = [1; xpi; xpi(1)^2; xpi(1) * xpi(2); xpi(2)^2];
        IPSPrimeX = [zeros(1, 1); ones(1, 1); zeros(1, 1); 2 * xpi(1); xpi(2); zeros(1, 1)];
        ISPPrimeY = [zeros(1, 1); zeros(1, 1); ones(1, 1); zeros(1, 1); xpi(1); 2 * xpi(2)];
        
        PhiIPS = P \ IPS;
        PhiDxIPS = P \ IPSPrimeX;
        PhiDyIPS = P \ ISPPrimeY;
        
        %bgFlowIPS = bgFlow(nodes);
        for k = 1:1
            Ux(j,k) = velocities{1,k,1}(nodes).' *  PhiIPS;
            Uy(j,k) = velocities{1,k,2}(nodes).' * PhiIPS;
        end
        %gammax(j) = bgFlowIPS * PhiDxIPS;
        %gammay(j) = bgFlowIPS * PhiDyIPS;
        xvec = 1/501:1/501:(500/501);
    end
    
    
end


paramSet = {[0.01  0.00000025  200]};
%xp = [2*ones(size(0.0:0.05:0.45)); 0.0:0.05:0.45];
xp = [2  ;0.25 ];
velocities = cell(length(paramSet), 1,2);

sol=@(x) zeros(size(x(:,1)));
f=@(x) -ones(size(x(:,1)));


fb=@(p) sol(p);

for i = 1:length(paramSet)
    parami = paramSet{i};
    
    fd=@(p) drectangle(p,0,4,0,1);
    fh=@(p) min(max(parami(1), parami(1)+1.0*drectangle(p, 1.99, 2.01, -0.1, 1.9).^3),0.1);
    [p,t]=distmesh2d(fd,fh,parami(1),[0, 0;4, 1],[0  0; 0  1;  4  1;  4   0]);
    maxWaveNum = parami(3);
    L = 4;
    % z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
    % z = circshift(z, -maxWaveNum/2+1);
    %waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;
    waveNumbers = -2*pi*(0:(maxWaveNum-1))/L;
    tic
    [bgFlow,~,~] = poissonSolver(p,t,f,fb);
    [u,p,t]=velocitySolveDiscDecompDiscTest(p,t, parami(2), waveNumbers, xp, bgFlow', L,2);
    for j = 1:size(xp,2)
        velocities(i,j,1) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,1)*2};
        velocities(i,j,2) = {(2*pi)^0.5*maxWaveNum/L*u(:,j,2)*2};
    end
    toc
    Uxd = zeros(500,1);
    Uyd = zeros(500,1);
    for j = 1:500
        xpi = [2; 1*j/501];
        xpn = whatTriangleIsThisPointIn(p,t,xpi);
        nodes = t(xpn, :);
        particleNodes{j} = nodes;
        % 6 by 6 matrix with rows: [ones; x; y; x^2; xy; y^2]:
        P = [ones(1, 6);
            p(:, nodes);
            p(1, nodes).^2;
            p(1, nodes) .* p(2, nodes);
            p(2, nodes).^2];
        IPS = [1; xpi; xpi(1)^2; xpi(1) * xpi(2); xpi(2)^2];
        IPSPrimeX = [zeros(1, 1); ones(1, 1); zeros(1, 1); 2 * xpi(1); xpi(2); zeros(1, 1)];
        ISPPrimeY = [zeros(1, 1); zeros(1, 1); ones(1, 1); zeros(1, 1); xpi(1); 2 * xpi(2)];
        
        PhiIPS = P \ IPS;
        PhiDxIPS = P \ IPSPrimeX;
        PhiDyIPS = P \ ISPPrimeY;
        
        %bgFlowIPS = bgFlow(nodes);
        for k = 1:1
            Uxd(j,k) = velocities{1,k,1}(nodes).' *  PhiIPS;
            Uyd(j,k) = velocities{1,k,2}(nodes).' * PhiIPS;
        end
        %gammax(j) = bgFlowIPS * PhiDxIPS;
        %gammay(j) = bgFlowIPS * PhiDyIPS;
        xvec = 1/501:1/501:(500/501);
    end
    
   figure
subplot(1,2,1)
hold on
plot(xvec-0.25,Uy,'lineWidth',2.0);
scatter(0 ,Uy(125) ,44 ,'r','filled')
%title('Line Slice Of Disturbance Velocity U-U_{str}','FontWeight','bold')
xlabel('y','FontSize',16,'FontWeight','bold')
ylabel('(U_1^{(1)} - U_{str})_y','FontSize',18,'FontWeight','bold')
%text(xd, yd, nl, 'FontSize',18, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
ax = gca;
ax.FontSize = 12;
subplot(1,2,2)
hold on
plot(xvec-0.25,Uyd,'lineWidth',2.0);
scatter(0 ,Uyd(125) ,44 ,'r','filled')
%title('Line Slice Of Corrected Disturbance Velocity V','FontWeight','bold')
xlabel('y','FontSize',16,'FontWeight','bold')
ylabel('(U_1^{(1)} - U_{str} - U_D)_y','FontSize',18,'FontWeight','bold')
%text(xd, yd, nl, 'FontSize',18, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
ax = gca;
ax.FontSize = 12;
end