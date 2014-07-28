function val=integral_base_a(a,b,c1,c2,nu,sc,n)
f=@(x)sc./(x.^2+nu^2).*(x-c1).*(x-c2);
[x,w]=lgwt(n,a,b);
val=w'*f(x);
end
