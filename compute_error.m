function [err1, err2]=compute_error(e1,e2,e1c,e2c,x,xc)
tic;
err1=compute_L2_error(e1,e1c,x,xc);
err2=compute_L2_error(e2,e2c,x,xc);
t=toc();
display (strcat(num2str(t)," for the computation of the total error"));

end

function [e1]=compute_L2_error(e, exact, x, xexact)



fexact=@(t)bsxfun(@times, exact(1:1:end-1), ...
	bsxfun(@times, bsxfun(@minus, xexact'(2:1:end), t), 1./(xexact'(2:1:end)-xexact'(1:1:end-1))));
qexact=@(t)bsxfun(@times, exact(2:1:end), ...
	bsxfun(@times, bsxfun(@minus, t, xexact'(1:1:end-1)), 1./(xexact'(2:1:end)-xexact'(1:1:end-1))));



quad_pts_current=@(a,b,quad_pts)bsxfun(@plus,bsxfun(@minus,b,a)./2*quad_pts,bsxfun(@plus,b,a)./2);
w_current=@(a,b,w)(b-a)./2*w;
[quad_pts, w]=lgwt(5,-1,1);
quad_pts=quad_pts';
w=w';
W=@(weights)w_current(xexact(1:1:end-1)',xexact(2:1:end)',weights);
weight_matrix=W(w);

Q=quad_pts_current(xexact(1:1:end-1)',xexact(2:1:end)',quad_pts);

FQ=zeros(size(Q));
c=1;
%hold on

for c=1:1:length(x)-1
	P=(xexact>=x(c))&(xexact<x(c+1));
	FQ(P,:)=e(c+1)*(Q(P,:)-x(c))./(x(c+1)-x(c))+e(c)*(x(c+1)-Q(P,:))./(x(c+1)-x(c));
end




F=bsxfun(@times, weight_matrix, (abs(FQ-fexact(Q)-qexact(Q))).^2);


e1=sqrt(sum(sum(abs(F))))




end
