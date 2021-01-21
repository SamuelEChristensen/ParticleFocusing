function Uhat=regStressletF(p,xp,gammax,gammay,epsilon, waveNum,L)
%pn is numberOfNodes.old
Z = -(L/2-L/waveNum):L/waveNum:L/2;
Z = circshift(Z, -waveNum/2+1);    
Uhat = fft(regStresslet(p, Z, xp, gammax, gammay, epsilon), waveNum,2);
end