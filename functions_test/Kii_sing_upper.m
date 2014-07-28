function val=Kii_sing_upper(x_l, x_c, x_r, alpha, nu)
hi=x_r-x_c;
val=0;
hiprev=x_c-x_l;
result=@(a,b,c)b^2./2-a^2./2-...
(2*c-1i*nu)*(b-a)+...
(c^2-nu^2-2i*c*nu)*1/2*log((b^2+nu^2)./(a^2+nu^2))+...
(2*nu*c+(c^2-nu^2)*1i)*(atan(b/nu)-atan(a/nu));
if (abs(alpha-x_l)<=1e-13)
val=0;
elseif alpha<=x_c
	val=1./(hiprev^2).*result(x_l,alpha,x_l);
elseif alpha>x_c
	val=1./(hiprev.^2)*result(x_l,x_c,x_l)+1./(hi^2).*result(x_c,alpha,x_r);


end
end
