function [e1, e2, A,B,D,x, e11, e21]=solve_mH1(dx,lambda, nu, uniform)
%mesh generation
L=20;
H=10;
if uniform
x=-L:dx:H;
else
phi=-pi/2:dx:pi/2;
 x=(H+L)*sin(phi)/2.0+(H-L)./2.0;
 %x=-L:dx/2.0:H/2;
 %x=[x (H./2+dx):dx:H];

end
tic
[A,B,D]=construct_block_matrix(x,nu,lambda);
time_passed=toc();
display (strcat(num2str(time_passed)," for the construction of the matrix"));
p=schur_diagonal_precond(A,B,B,D);
r=construct_rhs(x);



N=length(x);

tic
[e1, e2]=solve_block_system_gmres(A,B,D,r);%(A,B,B,D,r(1:N), r((N+1):end),p);
[e11,e21]=solve_block_system_naively(A,B,D,r(1:N), r((N+1):end),p);
time_passed=toc();
display(strcat(num2str(time_passed)," for the solution of the system of eqs"));

end




function a=alpha(x)

a=zeros(size(x));
%a(x<=-1)=1;
%P=(x<=3)&(x>-1);
%a(P)=-x(P);
%a(x>3)=-3;
%a=x;

a=x;

end

function delta=delta(x)

delta=zeros(size(x));
%delta(x<=-15)=0;
%P=(x>-15)&(x<=0);
%delta(P)=1+x(P)./15;
%P=(x>0)&(x<=3);
%delta(P)=1+x(P)./3;
%delta(x>3)=2;

delta=0;

end

%since the matrices are symmetric tridiagonal, each matrix constructor returns 
%a vector K of size $L\times 2$, where L is the size of the matrix
%the K(:,1) is the diagonal
%K(:,2) contains off-diagonal entries, with the last entry set to zero




%stiffness for P1
function K=K_lpl(x)
    L=length(x);
    K=zeros(L,2);

    K(1,1)=1/(x(2)-x(1));

    K(2:1:L-1,1)=1./(x(2:1:end-1)-x(1:1:end-2))+1./(x(3:1:end)-x(2:1:end-1));
    K(1:1:L-1,2)=-1./(x(2:1:end)-x(1:1:end-1));

    K(L,1)=1./(x(end)-x(end-1));
    

end

%mass matrix for P1
function M=mass_matrixP1(x)

    L=length(x);
    M=zeros(L,2);

    M(1,1)=(x(2)-x(1))/3;


    M(2:1:(L-1),1)=(x(2:1:end-1)-x(1:1:end-2))./3+(x(3:1:end)-x(2:1:end-1))./3;
    M(1:1:L-1,2)=(x(2:1:end)-x(1:1:end-1))./6;

    M(L,1)=(x(end)-x(end-1))./3;


end

%mass matrix for piecewise-constants
function M=mass_matrixL2(x)
    M=sparse(diagon(x(2:1:end)-x(1:1:end-1)));
end


%matrix M_{ij}=\int\limits_{-L}^{H}func_kernel(x)\phi_i(x) \phi_j(x)
function M=mass_matrixP1P1_scaled(x, func_kernel, nq)

    L=length(x);
    M=zeros(L,2);

    fleft_diagon=@(x_left,x_center,x)func_kernel(x).*(x-x_left).^2./(x_center-x_left).^2;
    fright_diagon=@(x_center,x_right,x)func_kernel(x).*(x-x_right).^2./(x_right-x_center).^2;


    fsubdiagon=@(x_left,x_right,x)func_kernel(x).*(x-x_left)./(x_right-x_left).*(x_right-x)./(x_right-x_left);


%nq is a number of quadrature points
    if(~nargin<3)
        nq=2;
    end

    %first value of the diagonal part is computed separately
    [quad_pts,w]=lgwt(nq(1),x(1),x(2));
    M(1,1)=w'*fright_diagon(x(1),x(2),quad_pts);



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
            M(j,1)=wc'*fleft_diagon(x(j-1),x(j),qpts);

            M(j-1,2)=wc'*fsubdiagon(x(j-1),x(j),qpts);
            %integrate 2nd part (for all values but the last one)
            if j<length(x)
                qpts=quad_pts_current(x(j),x(j+1),quad_pts);
                wc=w_current(x(j),x(j+1),w);
                M(j,1)=M(j,1)+wc'*fright_diagon(x(j),x(j+1),qpts);
            end
        end
       
    else
    %generate quadrature points for every subinterval separately
        for j=2:1:length(x)
            [qpts,wc]=lgwt(x(j-1),x(j), nq(j));
            M(j,1)=wc'*fleft_diagon(x(j-1),x(j),qpts);

            M(j-1,2)=wc'*fsubdiagon(x(j-1),x(j),qpts);
            %integrate 2nd part (for all values but the last one)
            if j<length(x)
                [qpts,wc]=lgwt(x(j+1),x(j),nq(j));
                M(j,1)=M(j,1)+wc'*fright_diagon(x(j),x(j+1),qpts);
            end
        end

    end
    
end

% a (presumably) faster implementaiton of mass_matrixP1P1_scaled
function M=mass_matrixP1P1_scaled_fast(x, func_kernel, nq)

L=length(x);

M=zeros(L,2);


fleft_diagon=@(x_left,x_center,x)bsxfun(@rdivide, func_kernel(x).*(bsxfun(@minus, x, x_left)).^2, (x_center-x_left).^2);
fright_diagon=@(x_center,x_right,x)bsxfun(@rdivide, func_kernel(x).*(bsxfun(@minus, x, x_right)).^2,(x_right-x_center).^2);


fsubdiagon=@(x_left,x_right,x)-func_kernel(x).*bsxfun(@rdivide,bsxfun(@minus,x,x_left).*bsxfun(@minus,x,x_right),(x_right-x_left).^2);

%nq is a number of quadrature points
if(nargin<3)
nq=2;
end
%diagon=zeros(1,length(x));
%subdiagon=zeros(1,length(x)-1);


%translating quadpoints
quad_pts_current=@(a,b,quad_pts)bsxfun(@plus,bsxfun(@minus,b,a)./2*quad_pts,bsxfun(@plus,b,a)./2);
w_current=@(a,b,w)(b-a)./2*w;


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
                 
M(1,1)=diagon_right(1);
M(2:1:L-1,1)=diagon_right(2:1:end)+diagon_left(1:1:end-1);
M(L,1)=diagon_left(end);
                 
func_interval_evaluator=@(q)fsubdiagon(x(1:1:end-1)',x(2:1:end)',q);

Mleft(1:1:(length(x)-1),1:end)=func_interval_evaluator(Q);
result=bsxfun(@times,weight_matrix,Mleft);
M(1:1:(L-1),2)=(sum(result.'))';
save('M.mat','M');


end
%returns the symmetrical tridiagonal matrices for the block matrix
%(A B)
%(C D)
              %A(:,1) are the $L$ diagonal entries of the matrix A
              %A(:,2) are the $L-1$ offdiagonal entries (and the last entry of the vector is set to zero)
%C is the complex conjugate of B
function [A,B,D]=construct_block_matrix(x,nu,lambda)
              L=length(x);
    
    Mdelta=mass_matrixP1P1_scaled_fast(x, @delta,4);

    Malpha=mass_matrixP1P1_scaled_fast(x, @alpha,4);
    
    
    M=mass_matrixP1(x);

    S=K_lpl(x);

                   A=Malpha+1i*nu*M;
                   B=1i*Mdelta;
                   D=S-(Malpha+1i*nu*M);
%adding a boundary condition
                   D(1,1)=D(1,1)-1i*lambda;

end

function b=construct_rhs(x)
    b=zeros(2*length(x),1);
    b(length(x)+1)=-2*1i*sqrt(2)*exp(1i*sqrt(2)*(-22));
end





%multiplies (A-BD^{-1}C)x

function r=multiply_schur_complement(A,B,C,D,x)
                    f1=multiply_tri_sym(C,x);
                    f2=D\f1;
                    f=multiply_tri_sym(B,f2);
                    f0=multiply_tri_sym(A,x);
                    r=f0-f;              
end


%computes the diagonal elements of the inverse of the tridiagonal matrix D
%(Rybicky-Hummer http://www.lanl.gov/DLDSTP/fast/diagonal.pdf)
function lambda=compute_diag_inv(D)
  b=D(:,1);
  n=length(D);
  a=-[0; D(1:1:n-1,2)];
  c=-[D(1:1:n-1,2); 0];
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
  p=A(:,1)-B(:,1).*dd.*C(:,1);
end
  
  
%solves the block system
%Ax+By=a
  %B'x+Dy=br
%by Schur complement and GMRES
function [x,y]=solve_block_system_naively(A,B,D,a,br)
  K=zeros(length(A)+length(B));
  L=length(A);
  K(1:1:L,1:1:L)=sparse(diag(A(:,1))+diag(A(1:1:L-1,2),-1)+diag(A(1:1:L-1,2),1));
  K((L+1):1:(2*L),1:1:L)=sparse(diag(B(:,1))+diag(B(1:1:L-1,2),-1)+diag(B(1:1:L-1,2),1));
  K(1:1:L,(L+1):1:(2*L))=sparse(diag(B(:,1))+diag(B(1:1:L-1,2),-1)+diag(B(1:1:L-1,2),1));
  K((L+1):1:(2*L), (L+1):1:(2*L))=sparse(diag(D(:,1))+diag(D(1:1:L-1,2),-1)+diag(D(1:1:L-1,2),1));
  r=[a;br];
  R=K\r;
  x=R(1:1:L);
  y=R((L+1):1:(2*L));
end


function [x,y]=solve_block_system_gmres(A,B,D,rhs)
mv=@(x)block_mv(A,B,B,D,x);
[h, flag, relres, iter, resvec] = gmres(mv, rhs, 5, 1e-5, 1000);
save("h.mat", 'h');
%for debugging purposes
display "gmres res: ";
display(flag);
display(iter);
display(relres);
Lx=length(A);
Ly=length(B);
x=h(1:1:Lx);
y=h((Lx+1):1:(Lx+Ly));
end



function [x,y]=solve_block_system(A,B,C,D,a,br, diagon_precond)
  
n=length(br);
%transform D into a format known to octave
Dmat=spdiags([D(:,2) D(:,1) circshift(D(:,2),1)], -1:1, length(D), length(D));

%first solve the system (A-BD^{-1}C)x=a-BD^{-1}br
%1. form the rhs
f=Dmat\br;
h=multiply_tri_sym(B,f);
rhs=a-h;

 
  
%2. solve the system with GMRES

mv=@(x)multiply_schur_complement(A,B,C,Dmat,x);
prec_mul=@(x)diagon_precond.*x;

%maxit num : 40
%gmres iteration restart: 5 is numerically good for the case without the preconditioner as well

[x, flag, relres, iter, resvec] = gmres(mv, rhs, 5, 1e-12, 40, prec_mul);



%next solve the remaining system Dy=-Cx+br;
h=multiply_tri_sym(C,x);
rhs=br-h;
y=Dmat\rhs;

end

function r=block_mv(A,B,C,D,h)
Lx=size(A);
Ly=size(C);
Lx=Lx(1);
Ly=Ly(1);
r1=multiply_tri_sym(A,h(1:1:Lx))+multiply_tri_sym(B, h((Lx+1):1:(Lx+Ly)));
r2=multiply_tri_sym(C,h(1:1:Lx))+multiply_tri_sym(D, h((Lx+1):1:(Lx+Ly)));
r=[r1;r2];
end


  

  
  

























