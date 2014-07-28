function val=Ksubdiag_x_upper(xi_l, xi_r, alpha)
val=0;
hi=xi_r-xi_l;

if abs(alpha-xi_l)<=1e-13
 val=0;
else if abs(alpha-xi_r)<=1e-13
 val=hi^2/12+xi_l*hi/6;
else
 p=(alpha-xi_l)/hi;
 val=hi^2*(p^3/3-p^4/4)+hi*xi_l*(p^2/2-p^3/3);
end



end
