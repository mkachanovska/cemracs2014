function [err1, err2, err2h1]=compute_error_exact(e1,e2,x,exact1, exact2, exact2der)
%L2 error
err1=compute_L2_error(e1,exact1,x);
err2=compute_L2_error(e2, exact2,x);
err2h1=compute_H1_error(e2, exact2der, x);
err2h1=sqrt(err2h1^2+err2^2);
end

function [e1]=compute_L2_error(e, exact, x)
e1=0;
f=@(t)bsxfun(@times, e(1:1:end-1), bsxfun(@times, bsxfun(@minus, x'(2:1:end), t), 1./(x'(2:1:end)-x'(1:1:end-1))));
q=@(t)bsxfun(@times, e(2:1:end), bsxfun(@times, bsxfun(@minus, t, x'(1:1:end-1)), 1./(x'(2:1:end)-x'(1:1:end-1))));


quad_pts_current=@(a,b,quad_pts)bsxfun(@plus,bsxfun(@minus,b,a)./2*quad_pts,bsxfun(@plus,b,a)./2);
w_current=@(a,b,w)(b-a)./2*w;
[quad_pts, w]=lgwt(5,-1,1);
quad_pts=quad_pts';
w=w';
W=@(weights)w_current(x(1:1:end-1)',x(2:1:end)',weights);
weight_matrix=W(w);

Q=quad_pts_current(x(1:1:end-1)',x(2:1:end)',quad_pts);

F=bsxfun(@times, weight_matrix, (f(Q)+q(Q)-exact(Q)).^2);


e1=sqrt(sum(sum(F)))


end

function [e2]=compute_H1_error(e, exact, x)
e2=0;
f=@(t)bsxfun(@times, e(1:1:end-1), -1./(x'(2:1:end)-x'(1:1:end-1)));
q=@(t)bsxfun(@times, e(2:1:end), 1./(x'(2:1:end)-x'(1:1:end-1)));


quad_pts_current=@(a,b,quad_pts)bsxfun(@plus,bsxfun(@minus,b,a)./2*quad_pts,bsxfun(@plus,b,a)./2);
w_current=@(a,b,w)(b-a)./2*w;
[quad_pts, w]=lgwt(6,-1,1);
quad_pts=quad_pts';
w=w';
W=@(weights)w_current(x(1:1:end-1)',x(2:1:end)',weights);
weight_matrix=W(w);

Q=quad_pts_current(x(1:1:end-1)',x(2:1:end)',quad_pts);
A= (f(Q)+q(Q)-exact(Q)).^2;

F=bsxfun(@times, weight_matrix, A);
P=F;


e2=sqrt(sum(sum(F)))


%B=zeros(size(A));
%e23=0;
%for k=1:1:length(x)-1
	%compute the error on the interval x(k) x(k+1)
	%the approximation given by P1 elements
%	f=@(t)e(k)*(-1)./(x(k+1)-x(k))+e(k+1)./(x(k+1)-x(k));
%	h=@(t)(f(t)-exact(t)).^2;
%		[qpts,w]=lgwt(6,x(k), x(k+1));
%	e23=e23+sum(w.*h(qpts));
%end
%e23=sqrt(e23)-e2
end

