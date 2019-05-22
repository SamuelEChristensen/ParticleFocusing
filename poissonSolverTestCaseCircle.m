%test case for oseenSolver on a circle.
 fd=@(p) sqrt(sum(p.^2,2))-1;
 N=@(x) exp(-1/2*(x(:,1).^2+x(:,2).^2));

initialLengths=[0.25,0.2,0.15,0.1,0.05];

error=zeros(length(initialLengths),2);


sol=@(x) [x(:,1).^2.*x(:,2)];
f=@(x) [-2*x(:,2)];

%sol=@(x,k) [zeros(size(x,1),1), ones(size(x,1),1),ones(size(x,1),1)];
%f=@(x,k) [x(:,1),-x(:,2).^2.*x(:,1),3*ones(size(x,1),1)];


fb=@(p) sol(p);

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[0,0]);
fullsol=sol(p);
    %profile on
    [Uwn,pold,told]=poissonSolver(p,t,f,fb);
    %profile viewer
    error(i,1)=norm(Uwn(1:length(pold))-fullsol)/length(p);
    error(i,2)=length(p);
    %plotFESol(p,t,Uwn(1:length(pold)))
    %plotFESol(p,t,abs(Uwn{2,2}-fullsol(:,2)))
    %plotStokesSol(p,t,Uwn);
    
end
figure
plot(log(error(:,2)),log(error(:,1)))