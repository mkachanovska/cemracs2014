function val=Kii_const_upper_a(xi_l, xi_c, xi_r, alpha)
val=0;
hi=xi_r-xi_c;
hiprev=xi_c-xi_l;

f=@(x)(1-(x-xi_c)/hi).^2;
g=@(x)((x-xi_l)/hiprev).^2;

if alpha<=xi_c
  [x,w]=lgwt(3,xi_l, alpha);
  val=w'*g(x);
else
  [x,w]=lgwt(3,xi_c,alpha);
  val=w'*f(x)
  [x,w]=lgwt(3,xi_l,xi_c);
  val=val+w'*g(x)
end
end
