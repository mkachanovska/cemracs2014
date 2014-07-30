function [e1, e2, A,B,C,D,precond, x, M, G, Morig, Gorig]=solve_mH1(dx,lambda, nu)
%mesh generation
L=2;
H=10;

x=-L:dx:H;

tic
[A,B,C,D]=construct_block_matrix(x,nu,lambda);
toc
display ("for the construction of the matrix");

r=construct_rhs(x);



N=length(x);

tic
precond=schur_diagonal_precond(A,B,C,D);
diagon_precond=@(x)precond.*x;
[e1, e2]=solve_block_system(A,B,C,D,r(1:N), r((N+1):end),diagon_precond);
toc
display(" for the solution of the system of eqs");
end


function a=alpha(x)
if length(x)==1
    if x<=-1 a=1;
    elseif x<=3 a=-x;
    else a=-3;
    end
else
a=zeros(size(x));
a(x<=-1)=1;
P=(x<=3)&(x>-1);
a(P)=-x(P);
a(x>3)=-3;
end

end

function delta=delta(x)
if length(x)==1
    if x<=-1 delta=0;
    elseif x<=3 delta=0.5*(x+1);
    else delta = 2;
    end
else
delta=zeros(size(x));
delta(x<=-1)=0;
P=(x<=3)&(x>-1);
delta(P)=0.5*(x(P)+1);
delta(x>3)=2;
end



end



%stiffness for P1
function K=K_lpl(x)
    L=length(x);
    diagon=zeros(L,1);
    diagon_und=zeros(L-1,1);
    diagon(1)=1/(x(2)-x(1));

    diagon(2:1:end-1)=1./(x(2:1:end-1)-x(1:1:end-2))+1./(x(3:1:end)-x(2:1:end-1));
    diagon_und(1:1:end)=-1./(x(2:1:end)-x(1:1:end-1));

    diagon(end)=1./(x(end)-x(end-1));

    K=sparse(diag(diagon)+diag(diagon_und,-1)+diag(diagon_und,1));
end

%mass matrix for P1
function K=mass_matrixP1(x)

    L=length(x);
    diagon=zeros(L,1);
    diagon_und=zeros(L-1,1);
    diagon(1)=(x(2)-x(1))/3;


    diagon(2:1:end-1)=(x(2:1:end-1)-x(1:1:end-2))./3+(x(3:1:end)-x(2:1:end-1))./3;
    diagon_und(1:1:end)=(x(2:1:end)-x(1:1:end-1))./6;

    diagon(end)=(x(end)-x(end-1))./3;


    K=sparse(diag(diagon)+diag(diagon_und, -1)+diag(diagon_und, 1));
end

%mass matrix for piecewise-constants
function M=mass_matrixL2(x)
    M=sparse(diagon(x(2:1:end)-x(1:1:end-1)));
end


%matrix M_{ij}=\int\limits_{-L}^{H}func_kernel(x)\phi_i(x) \phi_j(x)
function M=mass_matrixP1P1_scaled(x, func_kernel, nq)

    fleft_diagon=@(x_left,x_center,x)func_kernel(x).*(x-x_left).^2./(x_center-x_left).^2;
    fright_diagon=@(x_center,x_right,x)func_kernel(x).*(x-x_right).^2./(x_right-x_center).^2;


    fsubdiagon=@(x_left,x_right,x)func_kernel(x).*(x-x_left)./(x_right-x_left).*(x_right-x)./(x_right-x_left);


%nq is a number of quadrature points
    if(~nargin<3)
        nq=2;
    end
    diagon=zeros(length(x),1);
    subdiagon=zeros(length(x)-1,1);

    %first value of the diagonal part is computed separately
    [quad_pts,w]=lgwt(nq(1),x(1),x(2));
    diagon(1)=w'*fright_diagon(x(1),x(2),quad_pts);



    %translating quadpoints
    quad_pts_current=@(a,b,quad_pts)(b-a)./2*quad_pts+(a+b)./2;
    w_current=@(a,b,w)(b-a)./2*w;


    if length(nq)==1
        [quad_pts,w]=lgwt(nq, -1, 1);

        for j=2:1:length(x)
            %integrate 1st part of diagonal ('left hat') and the subdiagonal
            %1) translate the quadrature points
            qpts=quad_pts_current(x(j-1),x(j),quad_pts);
            wc=w_current(x(j-1),x(j),w);
            %2) integrate 
            diagon(j)=wc'*fleft_diagon(x(j-1),x(j),qpts);

            subdiagon(j-1)=wc'*fsubdiagon(x(j-1),x(j),qpts);
            %integrate 2nd part (for all values but the last one)
            if j<length(x)
                qpts=quad_pts_current(x(j),x(j+1),quad_pts);
                wc=w_current(x(j),x(j+1),w);
                diagon(j)=diagon(j)+wc'*fright_diagon(x(j),x(j+1),qpts);
            end
        end
       
    else
    %generate quadrature points for every subinterval separately
        for j=2:1:length(x)
            [qpts,wc]=lgwt(x(j-1),x(j), nq(j));
            diagon(j)=wc'*fleft_diagon(x(j-1),x(j),qpts);

            subdiagon(j-1)=wc'*fsubdiagon(x(j-1),x(j),qpts);
            %integrate 2nd part (for all values but the last one)
            if j<length(x)
                [qpts,wc]=lgwt(x(j+1),x(j),nq(j));
                diagon(j)=diagon(j)+wc'*fright_diagon(x(j),x(j+1),qpts);
            end
        end

    end
    
    M=sparse(diag(diagon)+diag(subdiagon,-1)+diag(subdiagon,1));
end

% a (presumably) faster implementaiton of mass_matrixP1P1_scaled
function M=mass_matrixP1P1_scaled_fast(x, func_kernel, nq)

fleft_diagon=@(x_left,x_center,x)bsxfun(@rdivide, func_kernel(x).*(bsxfun(@minus, x, x_left)).^2, (x_center-x_left).^2);
fright_diagon=@(x_center,x_right,x)bsxfun(@rdivide, func_kernel(x).*(bsxfun(@minus, x, x_right)).^2,(x_right-x_center).^2);


fsubdiagon=@(x_left,x_right,x)-func_kernel(x).*bsxfun(@rdivide,bsxfun(@minus,x,x_left).*bsxfun(@minus,x,x_right),(x_right-x_left).^2);


%nq is a number of quadrature points
if(~nargin<3)
nq=2;
end

diagon=zeros(1,length(x));
subdiagon=zeros(1,length(x)-1);


%translating quadpoints
quad_pts_current=@(a,b,quad_pts)bsxfun(@plus,bsxfun(@minus,b,a)./2*quad_pts,bsxfun(@plus,b,a)./2);
w_current=@(a,b,w)bsxfun(@minus,b,a)./2*w;


[quad_pts,w]=lgwt(nq, -1, 1);
quad_pts=quad_pts';
w=w';
W=@(weights)w_current(x(1:1:end-1)',x(2:1:end)',weights);

weight_matrix=W(w);

Q=quad_pts_current(x(1:1:end-1)',x(2:1:end)',quad_pts);



%form the matrix (length(x)-1)\times nq
%consisting of the evaluations of fleft_diagon in the quad point inside the intervals [x_i,x_i+1]
func_interval_evaluator=@(q)fleft_diagon(x(1:1:end-1)',x(2:1:end)',q);
Mleft=zeros(length(x)-1,nq);
Mleft(1:1:(length(x)-1),:)=func_interval_evaluator(Q);

result=bsxfun(@times,weight_matrix,Mleft);

diagon_left=sum(result.');

                
%similarly for fright_diagon
func_interval_evaluator=@(q)fright_diagon(x(1:1:end-1)',x(2:1:end)',q);
Mleft(1:1:(length(x)-1),:)=func_interval_evaluator(Q);
result=bsxfun(@times,weight_matrix,Mleft);
diagon_right=sum(result.');
                 
diagon(1)=diagon_right(1);
diagon(2:1:end-1)=diagon_right(2:1:end)+diagon_left(1:1:end-1);
diagon(end)=diagon_left(end);
                 
func_interval_evaluator=@(q)fsubdiagon(x(1:1:end-1)',x(2:1:end)',q);

Mleft(1:1:(length(x)-1),1:end)=func_interval_evaluator(Q);
result=bsxfun(@times,weight_matrix,Mleft);
subdiagon=sum(result.');

M=sparse(diag(diagon)+diag(subdiagon,-1)+diag(subdiagon,1));


end

function [A,B,C,D]=construct_block_matrix(x,nu,lambda)
    L=length(x);
    K=zeros(2*L,2*L);
    
    Mdelta=mass_matrixP1P1_scaled_fast(x, @delta,3);

    Malpha=mass_matrixP1P1_scaled_fast(x,@alpha,3);

    M=mass_matrixP1(x);

    S=K_lpl(x);

 %   K(1:L,1:L)=sparse(-Malpha-1i*nu*M);
 %   K(1:L,(L+1):(2*L))=sparse(-1i*Mdelta);
 %   K((L+1):(2*L),1:L)=sparse(1i*Mdelta);
 %   K((L+1):(2*L),(L+1):(2*L))=sparse(S-Malpha-1i*nu*M);
    %adding a boundary condition
 %   K(L+1,L+1)=K(L+1,L+1)-1i*lambda;

A=sparse(-Malpha-1i*nu*M);
B=sparse(-1i*Mdelta);
C=sparse(1i*Mdelta);
D=sparse(S-Malpha-1i*nu*M);
%adding a boundary condition

D(1,1)=D(1,1)-1i*lambda;

end

function b=construct_rhs(x)
    b=zeros(2*length(x),1);
b(length(x)+1)=-2*1i*sqrt(2)*exp(1i*sqrt(2)*(-22));

end





%multiplies (A-BD^{-1}C)x
function r=multiply_schur_complement(A,B,C,D,x)

r1=A*x;

r_a1=C*x;

r_a2=D\r_a1;

r2=B*r_a2;
r=r1-r2;
end


%computes the diagonal elements of the inverse of the tridiagonal matrix D
%(Rybicky-Hummer http://www.lanl.gov/DLDSTP/fast/diagonal.pdf)
function lambda=compute_diag_inv(D)
  b=diag(D);
  a=-[0; diag(D,-1)];
  c=-[diag(D,+1); 0];
  n=length(D);
  e=zeros(n+1,1);
  d=zeros(n+1,1);
  e(n+1)=0;
  d(1)=0;
  for k=n:-1:1
  e(k)=a(k)./(b(k)-c(k)*e(k+1));
  end
  for k=2:1:n+1
  d(k)=c(k-1)./(b(k-1)-a(k-1)*d(k-1));
  end
  lambda=zeros(n,1);
  lambda(1:1:end)=(1-d(2:1:end).*e(2:1:end)).^(-1).*(b(1:1:end)-a(1:1:end).*d(1:1:end-1)).^(-1);
end


%forms the diagonal preconditioner for A-BD^{-1}C
function p=schur_diagonal_precond(A,B,C,D)
  dd=compute_diag_inv(D);
  ad=diag(A);
  bd=diag(B);
  cd=diag(C);
  p=ad-bd.*dd.*cd;
end
  
  
%solves the block system
%Ax+By=alpha
%Cx+Dy=beta
%by Schur complement and GMRES
function [x,y]=solve_block_system(A,B,C,D,alpha,beta,diagon_precond)
%first solve the system (A-BD^{-1}C)x=alpha-BD^{-1}beta
%1. form the rhs
f=D\beta;
rhs=alpha-B*f;
%2. solve the system with GMRES

mv=@(x)multiply_schur_complement(A,B,C,D,x);

%maxit num : 40
%gmres iteration restart: 5 is numerically good for the case without the preconditioner as well
  [x, flag, relres, iter, resvec] = gmres(mv, rhs, 5, [], 40,diagon_precond);

%for debugging purposes
display "gmres res: ";
display(flag);
display(iter);
display(relres);


%next solve the remaining system Dy=-Cx+beta;

rhs=-C*x+beta;
y=D\rhs;

end

  
  %solves the same system with the help of stabilized biconjugate gradients
  function [x,y]=solve_block_system_bicg(A,B,C,D,alpha,beta,diagon_precond)
  %first solve the system (A-BD^{-1}C)x=alpha-BD^{-1}beta
  %1. form the rhs
  f=D\beta;
  rhs=alpha-B*f;
  %2. solve the system with GMRES
  
  mv=@(x)multiply_schur_complement(A,B,C,D,x);
  
  [x, flag, relres, iter, resvec] =bicgstab(mv,rhs,[],[],diagon_precond);

  display "bicgstab res: ";
  display(flag);
  display(iter);
  display(relres);
  
  
  %next solve the remaining system Dy=-Cx+beta;
  
  rhs=-C*x+beta;
  y=D\rhs;
  
  end
  
  
  

























