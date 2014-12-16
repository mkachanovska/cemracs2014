function [e1, e2,M,x, Kstiff, Mmass]=solve_mH1(dx,lambda, nu, uniform)

%mesh generation
L=5; %20!!!
H=19;%10!!!
if uniform
x=-L:dx:H;
else
phi=-pi/2:dx:pi/2;
 x=(H+L)*sin(phi)/2.0+(H-L)./2.0;


end
tic
[A,B,D]=construct_block_matrix(x,nu,lambda);
time_passed=toc();
display (strcat(num2str(time_passed),' for the construction of the matrix'));
r=construct_rhs(x);

tic
[M,rhs]=permute(A,B,D,r);
%save('M.mat', 'M');
sol=M\rhs;

norm(M*sol-rhs)

[e1, e2]=permute_solution_back(sol);

time_passed=toc();
display(strcat(num2str(time_passed),' for the solution of the system of eqs'));

%auxiliary matrices, namely the mass matrix and the stiffness matrix
Kst=K_lpl(x);
Mm=mass_matrixP1(x);
Kstiff=spdiags([Kst(:,2) Kst(:,1) circshift(Kst(:,2),1)], -1:1, length(Kst), length(Kst));
Mmass=spdiags([Mm(:,2) Mm(:,1) circshift(Mm(:,2),1)], -1:1, length(Mm), length(Mm));
end


function [M,rhs]=permute(A,B,D,rhs_old)
n1=size(A);
n2=size(B);
band=zeros(n1(1)+n2(1), 4);%7-diagonal symmetric system; the 'band' contains the diagonals and subdiagonals
%column permutation
%e_new_{j}=e_{1,(j+1)/2} if j is odd
%e_new_{j}=e_{2,j/2} if j is even
rhs=zeros(length(rhs_old),1);
rhs(2)=rhs_old(n1(1)+1);
%forming the new symmetric matrix, with the row permutation being the same as the column permutation
%main diagonal
band(1:2:(end-1),1)=A(:,1);
band(2:2:end, 1)=D(:,1);
%subdiagonal
band(1:2:(end-1),2)=B(:,1);
band(2:2:(end),2)=B(:,2);
%subsubiagonal
band(1:2:(end-1),3)=A(:,2);
band(2:2:end, 3)=D(:,2);
%subsubsubdiagonal
band(1:2:(end-1),4)=B(:,2);

M=spdiags([band(:,4) band(:,3) band(:,2) band(:,1) ...
  circshift(band(:,2),1) circshift(band(:,3),2) circshift(band(:,4),3)], -3:3, n1(1)+n2(1), n1(1)+n2(1));


end
function [e1, e2]=permute_solution_back(sol)
e1=sol(1:2:(end-1));
e2=sol(2:2:end);

end


function n=n_e(x)
n=zeros(size(x));
[omega, omega_c]=main_parameters();

n(x<-0.5)=0.25;
P=(x>=-0.5 & x<=9);
n(P)=(1+x(P))/2;
n(x>9)=5;
%if omega_c==0
%    d=x;
%else
%    d=sqrt((2*omega^2+x).^2-4*(x+omega^2)*(omega^2-omega_c^2));
%end
%n=0.5*(2*omega^2+x+d);%for omega=1

%n(x<-1)=0.5*0.75;
%P=(x>=-1)&(x<5);
%n(P)=0.75*(1+x(P)/2);
%P=(x>=5);
%n(P)=7/2*0.75;
end

function [omega, omega_c]=main_parameters()
omega=1;
omega_c=sqrt(0.5);
end


function a=alpha(x)
[omega, omega_c]=main_parameters();
n=n_e(x);
a=omega^2*(1-n/(omega^2-omega_c^2));

end

function delta=delta(x)
[omega, omega_c]=main_parameters();
n=n_e(x);
delta=omega*omega_c*n/(omega^2-omega_c^2);
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
    if(nargin<3)
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
nq=2
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
%save('M.mat','M');


end
%returns the symmetrical tridiagonal matrices for the block matrix
%(A B)
%(C D)
              %A(:,1) are the $L$ diagonal entries of the matrix A
              %A(:,2) are the $L-1$ offdiagonal entries (and the last entry of the vector is set to zero)
%C is the complex conjugate of B
function [A,B,D]=construct_block_matrix(x,nu,lambda)
              L=length(x);
   %default: 4 quad pts
    nq=4; 
    Mdelta=mass_matrixP1P1_scaled_fast(x, @delta,nq);

    Malpha=mass_matrixP1P1_scaled_fast(x, @alpha,nq);
    
    
    M=mass_matrixP1(x);

    S=K_lpl(x);

                   A=Malpha+1i*nu*M;
                   B=1i*Mdelta;
                   D=S-(Malpha+1i*nu*M);
%adding a boundary condition
                   D(1,1)=D(1,1)-1i*lambda;

end

function b=construct_rhs(x)
%[omega, omega_c]=main_parameters();
    b=zeros(2*length(x),1);
    b(length(x)+1)=0.11;
%    b(length(x)+1)=-2*1i*sqrt(2)*exp(1i*sqrt(2)*(-22));
 %    b(length(x)+1)=(-airy(-20)*i*2-airy(1,-20));
 %   b(length(x)+1)=-i*2*0.2-0.1;
end




  

  
  

























