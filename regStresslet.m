function U=regStresslet(p,z, xp ,gammax,gammay,epsilon)
x = p(:,1)-xp(1);
y = p(:,2)-xp(2);
gxy = (gammax*x+gammay*y);
z2 = z.^2;
Z = [2*x.*z.*gxy + epsilon^2*gammax*z;
    2*y.*z.*gxy + epsilon^2*gammay*z; 
    2*z2.*gxy + epsilon^2*gxy];
U = -repmat(5./(4*((epsilon^2+(x.^2+y.^2))+z2).^(5/2)),3,1).*Z;
end