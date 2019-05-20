pv=[0,0;1,0;1,1;0.5,0.5;0,1;0,0];


initialLengths=[0.04];

error=zeros(length(initialLengths),2);
count=1;
N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));
sol=@(x,k) [x(:,1).*x(:,2).*N(x,k), ...
    -(x(:,1).^2-1).*N(x,k), ...
    ones(size(x,1),1)];
 %f=@(x,k) [x(:,1).*x(:,2).*(x(:,1).^2+x(:,1).^2-6).*N(x,k)-x(:,1).*N(x,k), ...
  % N(x,k).*(x(:,1).^4+x(:,1).^2.*(x(:,2)-7)-x(:,2).^2+4)-x(:,2).*N(x,k), ...
  % ones(size(x,1),1)+1i*(2*pi)^0.5*k*N(x,k)]-2*pi*k^2*sol(x,k);

%sol=@(x,k) [0*ones(size(x,1),1),0*ones(size(x,1),1),0*ones(size(x,1),1)];
%f=@(x,k) [1*ones(size(x,1),1),2*ones(size(x,1),1),3*ones(size(x,1),1)];

 
sol=@(x,k) 0*[x(:,1).^2.*x(:,2), -x(:,2).^2.*x(:,1), zeros(size(x,1),1)]*N([0,0],k);
f =@(x,k) [2*x(:,2), -2*x(:,1), zeros(size(x,1))]-2*pi*k^2*sol(x,k);

fb=@(x,k) sol(x,k);

for i=1:length(initialLengths)
    meshSize=initialLengths(i);
    fh=@(p) min(0.01+1*abs(dcircle(p,0.5,0.5,0)).^2,0.1);
    fd=@(p) dpoly(p,pv);
    [p,t]=distmesh2d(fd,fh,0.01,[0,0;1.5,1.5],pv);

    [Uwn,pold,told]=stokesSolver(p,t,f,fb);
    error(i,1)=norm([Uwn{4,1}, Uwn{4,2}, Uwn{4,3}]-sol(p,0))/length(p);
    error(i,2)=length(p);
    plotStokesSol(p,t,Uwn);
    
end