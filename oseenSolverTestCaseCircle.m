%test case for oseenSolver on a circle.
 fd=@(p) sqrt(sum(p.^2,2))-0.64;
%fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.6));
 N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));

%initialLengths=[0.2,0.15,0.1,0.075,0.05];
initialLengths=[0.03];
error=zeros(length(initialLengths),2);


%sol=@(x,k) [1*(1-exp(3*x(:,2).^2+2*x(:,1).^2)), ...
%    (x(:,1).^2-1).*N(x,k), ...
%    3*sin(x(:,1)-xp(1))+3*cos(x(:,2)-xp(2))];
% f=@(x,k) [x(:,1).*x(:,2).*(x(:,1).^2+x(:,1).^2-6).*N(x,k)-x(:,1).*N(x,k), ...
%    N(x,k).*(x(:,1).^4+x(:,1).^2.*(x(:,2)-7)-x(:,2).^2+4)-x(:,2).*N(x,k), ...
% +1i*(2*pi)^0.5*k*N(x,k)];
sig = 100;
sol=@(x,k) [0*x(:,2), 0*x(:,1), 0*x(:,1)+0*k];
f=@(x,k) [sig*x(:,2)+x(:,1), ...
    -2*x(:,1)-x(:,1).^3, x(:,1)];

%sol=@(x,k) [zeros(size(x,1),1), ones(size(x,1),1),ones(size(x,1),1)];
%f=@(x,k) [x(:,1),-x(:,2).^2.*x(:,1),3*ones(size(x,1),1)];


fb=@(x,k) sol(x,k);

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[0,0]);
fullsol=sol(p,0);
    %profile on
    [Uwn,pold,told]=oseenSolver(p,t,f,fb);
    %profile viewer
    error(i,1)=norm([Uwn{2,1}, Uwn{2,2}, Uwn{2,3}]-fullsol)/length(p);
    error(i,2)=length(p);
    %plotFESol(p,t,Uwn{2,1}-fullsol(:,1))
    %plotFESol(p,t,abs(Uwn{2,2}-fullsol(:,2)))
    plotStokesSol(p,t,Uwn);
    
end
figure
plot(log(error(:,2)),log(error(:,1)))