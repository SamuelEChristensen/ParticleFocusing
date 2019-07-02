%test case for oseenSolver with gaussian approximation of delta function
fd=@(p) sqrt(sum(p.^2,2))-1;% +1/10*sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
%fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.6));
%fd=@(p) drectangle(p,-1,1,-1,1);

xp = [0.1, 0];
sig = 0.05;
N = @(x,k) 1/((sig^2*pi*2)^0.5)*exp(-1/2*(((x(:,1)-xp(1)).^2+(x(:,2)-xp(2)).^2)/sig^2+sig^2*k^2));
%initialLengths=[0.2,0.15,0.1,0.075,0.05];
initialLengths=[0.05];
error=zeros(length(initialLengths),2);

sol=@(x,k) [0*x(:,2), 0*x(:,1), 0*x(:,1)+0*k];
gammax = xp(1);
gammay = xp(2);
f=@(x,k) [gammax*sig^2*pi*2*1i*k.*N(x,k),...
          gammay*sig^2*pi*2*1i*k.*N(x,k), -sig^2*pi*2*(gammax*(x(:,1)-xp(1)).*N(x,k)+gammay*(x(:,2)-xp(2)).*N(x,k))];

fb=@(x,k) sol(x,k);
maxWaveNum = 128;
L = 10; 
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1.5,-1.5;1.5,1],[-1,-1;-1,1 ; 1, 1; 1, -1]);
fullsol=sol(p,0);
    %profile on
    [Uwn,pold,told]=oseenSolver(p,t,f,fb, waveNumbers);
    U = ifftUwn(Uwn);
    %profile viewer
    
    %plotFESol(p,t,Uwn{2,1}-fullsol(:,1))
    %plotFESol(p,t,abs(Uwn{2,2}-fullsol(:,2)))
    plotStokesSol(p,t,U);
    
end
%figure
%plot(log(error(:,2)),log(error(:,1)))