function [t,y] = EM(f,T,ic,dt)

t = T(1):dt:T(2);
y = zeros(length(t),length(ic));
y(1,:) = ic;
for i = 2:length(t)
y(i,:) = y(i-1,:) + dt*f(t(i-1),y(i-1,:)')';
end