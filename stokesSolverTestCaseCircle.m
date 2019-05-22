%test case for stokesSolver on a rectangle.
  fd=@(p) sqrt(sum(p.^2,2))-1;%+0.1*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));
  xp=[0.1, 0.1];  
  N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));

initialLengths=[0.1];

error=zeros(length(initialLengths),2);
count=1;

sol=@(x,k) [zeros(size(x,1),1), ones(size(x,1),1),ones(size(x,1),1)];
f=@(x,k) [x(:,1),-x(:,2).^2.*x(:,1),3*ones(size(x,1),1)];


fb=@(x,k) sol(x,k);

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
[p,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[0,0]);
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