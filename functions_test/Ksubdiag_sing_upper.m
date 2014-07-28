function val=Ksubdiag_sing_upper(x_l,x_r,alpha,nu)
h=x_r-x_l;

atanval=atan(alpha/nu)-atan(x_l/nu);
logval=log((alpha^2+nu^2)./(x_l^2+nu^2))
val=-1./h^2*(1i*(nu*(alpha-x_l)-(nu^2-x_r*x_l)*atanval-(x_r+x_l)./2*nu*logval)+...
             (alpha^2/2-x_l^2/2)+(x_r*x_l-nu^2)/2*logval-(x_l+x_r)*(alpha-x_l)+nu*(x_r+x_l)*atanval);

end
