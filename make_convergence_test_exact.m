function [e1l2, e2l2, e2h1]=make_convergence_test_exact
h=2e-6*2.^[0:2:20];
e1l2=zeros(length(h),1 );
e2l2=zeros(length(h),1 );
e2h1=zeros(length(h),1 );

f=@(x)airy(x);
g=@(x)airy(1,x);
ex=@(x)airy(x).*(-1i).*sqrt(x.^2+1).*sqrt(x.^2+1+x)./(x.^2+1);
for k=1:1:length(h)
 [e1, e2,M,x, Kstiff, Mmass]=solve_mH1_airy(h(end-k+1),2, 0,1);
 [e1l2(end-k+1), e2l2(end-k+1), e2h1(end-k+1)]=compute_error_exact(e1,e2,x,ex,f,g);
end
save('e1l2.mat', 'e1l2');
save('e2l2.mat', 'e2l2');
save('e2h1.mat', 'e2h1');
end
