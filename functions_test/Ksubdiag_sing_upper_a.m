function val=Ksubdiag_sing_upper_a(xi_l, xi_r, alpha, nu, n)
val=0;
hi=xi_r-xi_l;
[x,w]=lgwt(n,xi_l,alpha);
g=@(x)x./(x.^2+nu^2).*(x-xi_l)/hi.*(1-(x-xi_l)/hi);
val=w'*g(x);
gnu=@(x)nu./(x.^2+nu^2).*(x-xi_l)/hi.*(1-(x-xi_l)/hi);
val=val+1i*w'*gnu(x);
end
