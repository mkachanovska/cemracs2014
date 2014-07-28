function val=Kii_x_upper(xi_l, xi_c, xi_r, alpha)
val=0;
hi=xi_r-xi_c;
hiprev=xi_c-xi_l;
if abs(alpha-xi_l)<=1e-13
 val=0;
elseif abs(alpha-xi_c)<=1e-13
 val=hiprev^2/4+hiprev*xi_l/3;
elseif abs(alpha-xi_r)<=1e-13
 val=hiprev^2/4+hiprev*xi_l/3+hi^2*1/12+xi_c*hi/3;
elseif alpha<xi_c
  val=1/4*(alpha-xi_l)^4/hiprev^2+xi_l*1/3*(alpha-xi_l)^3/hiprev^2;
elseif alpha>xi_c
  p=(alpha-xi_c)/hi;
  val=hiprev^2/4+hiprev*1/3*xi_l+...
  (hi^2+xi_c*hi)*(1-(1-p)^3)/3-hi^2*(1-(1-p)^4)/4;;
end
end
