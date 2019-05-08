function Uhat=regStressletF(p,xp,gammax,gammay,epsilon,a)
L=max(sum(p.^2,2));
waveNum=64;
z = -(L/2-L/waveNum):L/waveNum:L/2;
U = regStresslet(p, z, xp, gammax, gammay, epsilon, a);
Uhat = fft(U, waveNum, 2);
end