function val=Kii_const_upper(xi_l, xi_c, xi_r, alpha)
val=0;
hi=xi_r-xi_c;
hiprev=xi_c-xi_l;
if (abs(alpha-xi_r)<=1e-13)
val=hi/3+hiprev/3;
elseif (abs(alpha-xi_l)<=1e-13)
val=0;
elseif (alpha<=xi_c)
val=(alpha-xi_l)^3/(3*hiprev^2);
elseif (alpha>xi_c)
val=hiprev/3+hi/3*(1-(1-(alpha-xi_c)/hi)^3);
end
end
