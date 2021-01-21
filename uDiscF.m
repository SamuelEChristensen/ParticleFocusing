function udf = uDiscF(p,x0,y0,gx,gy,waveNum,L)
%pn is numberOfNodes.old
Z = -(L/2-L/waveNum):L/waveNum:L/2;
Z = circshift(Z, -waveNum/2+1);    
udf =  fft(outerSolDisconinuity(p,Z',x0,y0,0,gx,gy),waveNum,2);

end