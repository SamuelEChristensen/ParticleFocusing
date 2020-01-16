function udf = uDiscF(p,x0,y0,gx,gy,waveNum)
%pn is numberOfNodes.old
L =4;
Z = -(L/2-L/waveNum):L/waveNum:L/2;
Z = circshift(Z, -waveNum/2+1);    
udf =  fft(outerSolDisconinuity(p,Z',x0,y0,0,gx,gy),waveNum,2);

end