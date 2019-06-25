function U = ifftUwn(Uwn)

maxWaveNum = length(Uwn);
U1 = zeros(length(Uwn{1,1}), maxWaveNum);
U2 = zeros(length(Uwn{1,1}), maxWaveNum);
U3 = zeros(length(Uwn{1,1}), maxWaveNum);
P = zeros(length(Uwn{1,4}), maxWaveNum);
U = cell(maxWaveNum,4);
for i = 1:maxWaveNum
    U1(:,i) = Uwn{i,1};
    U2(:,i) = Uwn{i,2};
    U3(:,i) = Uwn{i,3};
    P(:,i) = Uwn{i,4};
end
%U1 = circshift(U1, -maxWaveNum/2+1,2);
%U2 = circshift(U2, -maxWaveNum/2+1,2);
%U3 = circshift(U3, -maxWaveNum/2+1,2);
%P = circshift(P, -maxWaveNum/2+1,2);
U1 = ifft(U1,maxWaveNum,2);
U2 = ifft(U2,maxWaveNum,2);
U3 = ifft(U3,maxWaveNum,2);
P = ifft(P,maxWaveNum,2);

for i = 1:maxWaveNum  
    U(i,1) = {U1(:,i)};
    U(i,2) = {U2(:,i)};
    U(i,3) = {U3(:,i)};
    U(i,4) = {P(:,i)};
end