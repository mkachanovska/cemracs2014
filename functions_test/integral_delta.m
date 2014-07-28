function val=integral_delta(a,b,c1,c2,nu)
val=0;
alpha=1;
[x,w]=lgwt(4,a,b);
polynomials=@(x)2*alpha*(x-c1).*(x-c2)+(x+1i*nu).*(x-c1).*(x-c2);
val=w'*polynomials(x);

valx=integral_base_x(a,b,c1,c2,nu,1.0);
valc=integral_base(a,b,c1,c2,nu,nu);
val=val+(alpha^2-nu^2)*valx+1i*(alpha^2-nu^2)*valc;
val=val+2*alpha*1i*nu*valx-2*alpha*nu*valc;
end
