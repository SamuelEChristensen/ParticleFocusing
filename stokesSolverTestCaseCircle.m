%test case for stokesSolver on a rectangle.
  fd=@(p) sqrt(sum(p.^2,2))-1+0.1*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));
  xp=[0.1, 0.1];  
  N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));

initialLengths=[0.1];

error=zeros(length(initialLengths),2);
count=1;
sol=@(x,k) [1*(1-exp(3*x(:,2).^2+2*x(:,1).^2)), ...
    (x(:,1).^2-1).*N(x,k), ...
    3*sin(x(:,1)-xp(1))+3*cos(x(:,2)-xp(2))];
 f=@(x,k) [x(:,1).*x(:,2).*(x(:,1).^2+x(:,1).^2-6).*N(x,k)-x(:,1).*N(x,k), ...
    N(x,k).*(x(:,1).^4+x(:,1).^2.*(x(:,2)-7)-x(:,2).^2+4)-x(:,2).*N(x,k)+(x(:,1)).*(x(:,2))./((x(:,1)-xp(1)).^2+(x(:,2)-xp(2)).^4+0.00001), ...
  1/((x(:,1)-xp(1)).^2+(x(:,2)-xp(2)).^2+0.005)+1i*(2*pi)^0.5*k*N(x,k)];

%sol=@(x,k) [0*ones(size(x,1),1),0*ones(size(x,1),1),0*ones(size(x,1),1)];
%f=@(x,k) [1*ones(size(x,1),1),2*ones(size(x,1),1),3*ones(size(x,1),1)];


fb=@(x,k) sol(x,k);

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
[p,t]=distmesh2d(fd,@huniform,0.05,[-1,-1;1,1],[0,0]);
%fullsol=sol(p,0);
    profile on
    [Uwn,pold,told]=stokesSolver(p,t,f,fb);
    profile viewer
%    error(i,1)=norm([Uwn{2,1}, Uwn{2,2}, Uwn{2,3}]-fullsol)/length(p);
    error(i,2)=length(p);
    %plotFESol(p,t,Uwn{2,1}-fullsol(:,1))
  %  plotFESol(p,t,abs(Uwn{2,2}-fullsol(:,2)))
    plotStokesSol(p,t,Uwn);
    
end
figure
plot(log(error(:,2)),log(error(:,1)))