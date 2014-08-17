h0=1.0;
k=0:2:16
nu=2.^(-k);
href=2e-6;
l=1:2:21;
h=href*2.^l;
result_l2e1=zeros(length(nu), length(h));
result_l2e2=zeros(length(nu), length(h));
result_h1e2=zeros(length(nu), length(h));
result_pwe1=zeros(length(nu), length(h));
result_pwe2=zeros(length(nu), length(h));

for k=1:1:length(nu)
	%compute the reference solution for nu(k)
	[e1ref,e2ref,M,xref]=solve_mH1(href,sqrt(10.0),nu(k),1);
	fname=strcat('e1.', num2str(2e-6), '.', num2str(nu(k)), '.mat');
                save(fname, 'e1ref');
  	fname=strcat('e2.', num2str(2e-6), '.', num2str(nu(k)), '.mat');
                save(fname, 'e2ref');
	for l=1:1:length(h)
		%compute the solution for the current h(l)
		[e1,e2,M,x,Kstiff, Mmass]=solve_mH1(h(l),sqrt(10.0),nu(k),1);
		%compute the reference solutions on the coarser mesh
                fname=strcat('e1.', num2str(h(l)), '.', num2str(nu(k)), '.mat');
		save(fname, 'e1');
	        fname=strcat('e2.', num2str(h(l)), '.', num2str(nu(k)), '.mat');                save(fname, 'e2');
		
	%	[result_l2e1(k,l), result_l2e2(k,l)]=compute_error( e1, e2, e1ref, e2ref, x, xref);
	end
	save('result_l2e1.mat', 'result_l2e1');
	save('result_l2e2.mat', 'result_l2e2');
end

save('result_l2e1.mat', 'result_l2e1');
save('result_h1e2.mat', 'result_h1e2');
save('result_l2e2.mat', 'result_l2e2');
save('result_pwe1.mat', 'result_pwe1');
save('result_pwe2.mat', 'result_pwe2');






