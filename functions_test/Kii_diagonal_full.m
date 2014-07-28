function Kii_diagonal_full(x, i, nu)

supp=[x(i-1) x(i) x(i+1)];
val=0;
hi=supp(3)-supp(2);
hiprev=supp(2)-supp(1);
%simple cases
if supp(3)<=-1
	val=-(1+1i*nu)*(hiprev/3+hi/3);
elseif supp(1)>=3
	val=(4/(3+1i*nu)+3-1i*nu)*(hiprev/3+hi/3);
elseif supp(3)<=3 && supp(1)>=-1
	val=Kii_delta_alpha_diagon(supp(1),supp(2),supp(3), nu, supp(3));
elseif supp(1)<-1 && supp(2)>=-1
	valconst=Kii_const_upper(supp(1),supp(2),supp(3),-1);
	valrest=Kii_delta_alpha_diagon(supp(1),supp(2),supp(3),nu,supp(3)-...
		Kii_delta_alpha_diagon(supp(1),supp(2),supp(3),nu,-1);
	val=valconst+valrest;
elseif supp(3)>3 && supp(2)<=3
	valdelta=Kii_delta_alpha_diagon(supp(1),supp(2),supp(3),nu,3)+...
		 (4/(3+1i*nu)+3-1i*nu)*Kii_const_lower(supp(1),supp(2),supp(3),3);
end

end
