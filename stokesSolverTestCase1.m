%test case for stokesSolver on a rectangle.
%fd=@(p) drectangle(p,0,1,0,1);
N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));
fd=@(p) sqrt(sum(p.^2,2))-1;

initialLengths=[0.25, 0.2, 0.15, 0.1, 0.075, 0.05];
initialLengths=[0.1];

error=zeros(length(initialLengths),2);
count=1;

realSol=@(x,z) [0*N(x,z), ...
   -z*N(x,z), ...
  x(:,2).*N(x,z)];
sol =@(x,k) [zeros(size(x,1),1), ...
    -1i*k*N(x,k), ...
    x(:,2).*N(x,k)];
f = @(x,k) [zeros(size(x,1),1), ...
            -1i*k*N(x,k).*(x(:,1).^2 + x(:,2).^2 - 2) + k^2*1i*k*N(x,k), ...
            (x(:,2).^2+x(:,1).^2-4).*x(:,2).*N(x,k) - k^2*x(:,2).*N(x,k)];
            

%sol = @(x,k) [x(:,1).^2.*x(:,2),-x(:,1).*x(:,2).^2,0*ones(size(x,1),1)];
%f = @(x,k) [2*x(:,2),-2*x(:,1),0*ones(size(x,1),1)];


%sol = @(x,k) [1i*x(:,1)*k, 1i*x(:,2)*k, 2*ones(size(x,1),1)];
%f = @(x,k) [-1i*k^3*x(:,1), -1i*k^3*x(:,2), -2*k^2];

%sol=@(x,k) [x(:,1).^2.*(x(:,2)), -1*x(:,2).^2.*(x(:,1)), zeros(size(x,1),1)]+0*N([0,0],k);
%f =@(x,k) [2*x(:,2).*N(x,k)+x(:,2), -2*x(:,1)+x(:,1), zeros(size(x,1))]-2*pi*k^2*sol(x,k);
%sol=@(x,k) [sin(pi*x(:,1))+sin(pi*x(:,2)),-sin(pi*x(:,1))-sin(pi*x(:,2)), zeros(size(x,1),1)]*N([0,0],k);
fb=@(x,k) sol(x,k);

for i=1:length(initialLengths)
    meshSize=initialLengths(i);
     fh=@(p) min(0.01+0.3*abs(fd(p)),0.1);
     %[p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0; 0, 1; 1, 1; 1, 0]);

    %[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[0,0;1,1],[0,0; 0, 1; 1, 1; 1, 0]);
    [p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[0,0]);
    [Uwn,pold,told]=stokesSolver(p,t,f,fb);
    U = ifftUwn(Uwn);
    zwn = length(U)/2;
    for j = 1:length(U)
        error(i,1)=error(i,1)+norm([U{j,1}, U{j,2}, U{j,3}]-realSol(p,j-zwn))/length(p);
    end
    error(i,2)=length(p);
    %plotStokesSol(p,t,Uwn, 1);
    
end
figure
plot(log(error(:,2)),log(error(:,1)))