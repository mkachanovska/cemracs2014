function val=Ksubdiag_x_upper_a(xi_l, xi_r, alpha)
val=0;
hi=xi_r-xi_l;
[x,w]=lgwt(3,xi_l,alpha);
g=@(x)x.*(x-xi_l)/hi.*(1-(x-xi_l)/hi);
val=w'*g(x);
end
