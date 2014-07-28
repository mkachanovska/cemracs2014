function val=integral_deltaa_a(a,b,c1,c2,nu,n)
alpha=1;
f=@(x)(x+1i*nu).*(x+alpha).^2./(x.^2+nu^2).*(x-c1).*(x-c2);
[x,w]=lgwt(n,a,b);
val=w'*f(x);
end
