function udf = udzDiscF(p,x0,y0,gx,gy,waveNum,L, bgFlow)
%pn is numberOfNodes.old
Z = -(L/2-L/waveNum):L/waveNum:L/2;
Z = circshift(Z, -waveNum/2+1);  
ud = outerSolDisconinuity(p,Z',x0,y0,0,gx,gy);
udf =  fft(repmat(bgFlow-gx*p(:,1)-gy*p(:,2),3,1).*ud,waveNum,2);

end