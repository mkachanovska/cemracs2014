function val=integral_base(a,b,c1,c2,nu,sc)
%if abs(nu^2+a*b)>0
%  atanval=atan((b-a).*nu./(nu^2+a*b))
%else 
%  atanval=sign(b-a)*pi/2;
%end
atanval=atan(b/nu)-atan(a/nu);
logval=log((b^2+nu^2)./(a^2+nu^2))
val=sc*(b-a)-(nu^2-c1*c2)*sc/nu*atanval-(c1+c2)./2*sc*logval;
end
