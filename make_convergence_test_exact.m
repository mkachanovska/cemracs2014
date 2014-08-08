function [e1l2, e2l2, e2h1]=make_convergence_test_exact
p=0:1:7;
h=10.^(-7+p);
e1l2=zeros(length(h),1 );
e2l2=zeros(length(h),1 );
e2h1=zeros(length(h),1 );

f=@(x)airy(x);
g=@(x)airy(1,x);
ex=@(x)airy(x).*(-1i).*sqrt(x.^2+1).*sqrt(x.^2+1+x)./(x.^2+1);
for k=1:1:8
 [e1, e2,M,x, Kstiff, Mmass]=solve_mH1_airy(h(k),2, 0,1);
 [e1l2(k), e2l2(k), e2h1(k)]=compute_error_exact(e1,e2,x,ex,f,g);
end
save('e1l2.mat', 'e1l2');
save('e2l2.mat', 'e2l2');
save('e2h1.mat', 'e2h1');
end
