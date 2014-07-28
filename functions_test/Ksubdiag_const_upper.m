function val=Ksubdiag_const_upper(xi_l, xi_r, alpha)
val=0;
hi=xi_r-xi_l;
if abs(alpha-xi_r)<=1e-13
	val=hi/6;
elseif abs(alpha-xi_l)<=1e-13
	val=0;
else 
  val=hi^(-1)*((alpha-xi_l)^2/2-(alpha-xi_l)^3/(3*hi));
end
end
