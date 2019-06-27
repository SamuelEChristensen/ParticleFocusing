%test case for oseenSolver on a circle.
 fd=@(p) sqrt(sum(p.^2,2))-1;
%fd=@(p) ddiff(drectangle(p,0,0,-1,1),dcircle(p,0,0,0.6));
%fd=@(p) drectangle(p,0,1,0,1);
 N=@(x,k) exp(-1/2*(x(:,1).^2+x(:,2).^2+k^2));

initialLengths=[0.2, 0.15, 0.1, 0.075, 0.05];
%initialLengths=[0.09];
error=zeros(length(initialLengths),2);


realSol=@(x,z) [0*N(x,z), ...
   -z*N(x,z), ...
  x(:,2).*N(x,z)];
sol =@(x,k) [zeros(size(x,1),1), ...
    -1i*k*N(x,k), ...
    x(:,2).*N(x,k)];
% f = @(x,k) [zeros(size(x,1),1)...
%             + (1 + x(:, 1) + x(:, 2)) .* (1i * x(:, 2) .* -x(:,1) .* N(x,k)), ...
%             -1i*k*N(x,k).*(x(:,1).^2 + x(:,2).^2 - 2) + k^2*1i*k*N(x,k)...
%              + (1 + x(:, 1) + x(:, 2)) .* ( -x(:,2).^2 .* N(x,k)), ...
%             (x(:,2).^2+x(:,1).^2-4).*x(:,2).*N(x,k) - k^2*x(:,2).*N(x,k)...
%             - 1i*k*N(x,k) + (1 + x(:, 1) + x(:, 2)) .* (1i*k * x(:,2) .* N(x,k))];

f = @(x,k) [zeros(size(x,1),1)...
             + zeros(size(x,1),1), ...
            -1i*k*N(x,k).*(x(:,1).^2 + x(:,2).^2 - 2) + k^2*1i*k*N(x,k)...
              + (1 + x(:, 1) + x(:, 2)) .* (k^2*N(x,k)) , ...
            (x(:,2).^2+x(:,1).^2-4).*x(:,2).*N(x,k) - k^2*x(:,2).*N(x,k)...
            - 1i*k*N(x,k) + (1 + x(:, 1) + x(:, 2)) .* (1i*k * x(:,2) .* N(x,k))];

fb=@(x,k) sol(x,k);


maxWaveNum = 128;
L = 10; 
z = linspace(-L/2 + L/maxWaveNum, L/2, maxWaveNum);
z = circshift(z, -maxWaveNum/2+1);
waveNumbers = -2*pi*[0:maxWaveNum/2, -(maxWaveNum/2-1):-1]/L;

for i=1:length(initialLengths)
%fh=@(p) min(0.05+0.21*abs(dcircle(p,xp(1),xp(2),0)),0.1);
%[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[0,0]);
[p,t]=distmesh2d(fd,@huniform,initialLengths(i),[-1,-1;1,1],[0,0]);
    %profile on
    [Uwn,pold,told]=oseenSolver(p,t,f,fb, waveNumbers);
    %profile viewer
    U = ifftUwn(Uwn);
    zwn = length(Uwn)/2;
    %for j = 1:length(Uwn)
    %    error(i,1)=error(i,1)+norm([Uwn{j,1}, Uwn{j,2}, Uwn{j,3}]-sol(p,waveNumbers(j)))/length(p);
    %end
    for j = 1:length(Uwn)
        error(i,1)=error(i,1)+norm( (2*pi)^0.5*maxWaveNum/L*[U{j,1}, U{j,2}, U{j,3}]-realSol(p,z(j)))/length(p);
    end
    error(i,2)=length(p);
    
end
figure
plot(log(error(:,2)),log(error(:,1)))