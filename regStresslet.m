function U=regStresslet(p,z, xp ,gammax,gammay,epsilon,a)
x = p(:,1)-xp(1);
y = p(:,2)-xp(2);
Z = [2*x.*z.*(gammax*x+gammay*y) + epsilon^2*gammax*z,...
    2*y.*z.*(gammax*x+gammay*y) + epsilon^2*gammay*z, ...
    2*z.^2.*(gammax*x+gammay*y) + epsilon^2*(gammax*x+gammay*y)];
Z = reshape(Z,length(p),length(z),3);
U = -5*a^3./(4*(epsilon^2+x.^2+y.^2+z.^2).^(5/2)).*Z;
end