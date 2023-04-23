

X = rs(:,1:100);
Y = [X ];
nx = size(X,1)
%% augment matrix with mirror images to enforce symmetry/anti-symmetry
for k=1:size(X,2)
    %xflip = reshape(flipud(reshape(X(:,k),nx,ny)),nx*ny,1);
    %Y(:,k+size(X,2)) = xflip;
end

%% compute mean and subtract;
Uavg = mean(Y,2);
%f1 = plotCylinder(reshape(Uavg,nx,ny),nx,ny);  % plot average wake

%% compute POD after subtracting mean (i.e., do PCA)
[UPSI,US,UV] = svd(Y-Uavg*ones(1,size(Y,2)),'econ');
% PSI are POD modes
figure
semilogy(diag(US)./sum(diag(US))); % plot singular vals

for k=1:4  % plot first four POD modes
  %  f1 = plotCylinder(reshape(UPSI(:,k),nx,ny),nx,ny);
end