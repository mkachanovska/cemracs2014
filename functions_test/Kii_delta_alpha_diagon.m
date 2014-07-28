function val=Kii_delta_alpha_upperon(xi_l, xi_c, xi_r, nu, alpha)

hi=xi_r-xi_c;
hiprev=xi_c-xi_l;

f=@(xst, xend, c1, h)1/(4*h^2)*integral_delta(xst,xend,c1,c1,nu);
val=0;
if alpha<=xi_c
 val=f(xi_l,xi_c,xi_l,hiprev)+Kii_x_upper(xi_l, xi_c, xi_r, alpha)-1i*nu*Kii_const_upper(xi_l,xi_c,xi_r,alpha);
else 
 val=f(xi_l,xi_c,xi_l,hiprev)+f(xi_l,xi_c,xi_l,hi)+Kii_x_upper(xi_l, xi_c, xi_r, alpha)-1i*nu*Kii_const_upper(xi_l,xi_c,xi_r,alpha);
end
end
