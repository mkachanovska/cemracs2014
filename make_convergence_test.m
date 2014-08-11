h0=1.0;
k=0:1:15
nu=2.^(-k);
href=1e-6;
l=1:1:20;
h=href*2.^l;
result_l2e1=zeros(length(nu), length(h));
result_l2e2=zeros(length(nu), length(h));
result_h1e2=zeros(length(nu), length(h));
result_pwe1=zeros(length(nu), length(h));
result_pwe2=zeros(length(nu), length(h));

for k=1:1:16
	%compute the reference solution for nu(k)
	[e1ref,e2ref,M,xref]=solve_mH1(href,2.0,nu(k),1);
	for l=1:1:20
		%compute the solution for the current h(l)
		[e1,e2,M,x,Kstiff, Mmass]=solve_mH1(h(l),2.0,nu(k),1);
		%compute the reference solutions on the coarser mesh
	
		
		[result_l2e1(k,l), result_l2e2(k,l)]=compute_error( e1, e2, e1ef, e2ref, x, xref);
	end
end

save('result_l2e1.mat', 'result_l2e1');
save('result_h1e2.mat', 'result_h1e2');
save('result_l2e2.mat', 'result_l2e2');
save('result_pwe1.mat', 'result_pwe1');
save('result_pwe2.mat', 'result_pwe2');






