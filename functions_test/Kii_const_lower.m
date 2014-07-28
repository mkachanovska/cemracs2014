function val=Kii_const_lower(xi_l, xi_c, xi_r, alpha)
val=0;
hi=xi_r-xi_c;
hiprev=xi_c-xi_l;
val=hiprev/3+hi/3-Kii_const_upper(xi_l,xi_c,xi_r,alpha);
end

