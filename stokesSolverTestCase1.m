%test case for stokesSolver on a rectangle.
fd=@(p) drectangle(p,-1,1,-1,1);
N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));
%fd=@(p) sqrt(sum(p.^2,2))-1;

initialLengths=[0.3,0.25, 0.2, 0.15, 0.1];
%initialLengths=[0.075];

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
            
fb=@(x,k) sol(x,k);

profile clear
profile on

maxWaveNum = 128;
L = 10; 
%z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
%waveNumbers = -2*pi*(-(maxWaveNum/2-1):(maxWaveNum/2))/L;
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;



for i=1:length(initialLengths)
    meshSize=initialLengths(i);
     fh=@(p) min(0.01+0.3*abs(fd(p)),0.1);
    %[p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0; 0, 1; 1, 1; 1, 0]);

    %[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[0,0;1,1],[0,0; 0, 1; 1, 1; 1, 0]);
    [p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[1,1;-1,1;-1,-1;1,-1]);
    [Uwn,pold,told]=stokesSolver(p,t,f,fb, waveNumbers);
    U = ifftUwn(Uwn);
    zwn = length(Uwn)/2;
    %for j = 1:length(Uwn)
    %    error(i,1)=error(i,1)+norm([Uwn{j,1}, Uwn{j,2}, Uwn{j,3}]-sol(p,waveNumbers(j)))/length(p);
    %end
    for j = 1:length(Uwn)
        error(i,1)=error(i,1)+norm( (2*pi)^0.5*maxWaveNum/L*[U{j,1}, U{j,2}, U{j,3}]-realSol(p,z(j)))/length(p);
    end
    error(i,2)=length(p);
    %plotStokesSol(p,t,Uwn, 1);
    
end
profile viewer
figure
plot(log(error(:,2)),log(error(:,1)))