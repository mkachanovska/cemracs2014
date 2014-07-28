function val=Kii_sing_upper_a(xi_l, xi_c, xi_r, alpha, nu,n)
val=0;
hi=xi_r-xi_c;
hiprev=xi_c-xi_l;

f=@(x)x./(x.^2+nu^2).*(1-(x-xi_c)/hi).^2;
g=@(x)x./(x.^2+nu^2).*((x-xi_l)/hiprev).^2;

fnu=@(x)nu./(x.^2+nu.^2).*(1-(x-xi_c)/hi).^2;
gnu=@(x)nu./(x.^2+nu.^2).*((x-xi_l)/hiprev).^2;

if alpha<=xi_c
  [x,w]=lgwt(n,xi_l, alpha);
  val=w'*g(x)+1i*w'*gnu(x);
else
  [x,w]=lgwt(n,xi_c,alpha);
  val=w'*f(x)+1i*w'*fnu(x);
  [x,w]=lgwt(n,xi_l,xi_c);
  val=val+w'*g(x)+1i*w'*gnu(x);
end
end
