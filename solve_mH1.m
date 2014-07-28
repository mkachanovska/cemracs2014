function [u,K,x]=solve_mH1(dx,lambda, nu)
%mesh generation
L=2;
H=10;

x=-L:dx:H;

tic
K=construct_block_matrix(x,nu,lambda);
r=construct_rhs(x);
toc
disp( "formed the matrix and the RHS");
tic
u=K\r;
toc
end


function a=alpha(x)

    if x<=-1 a=1;
    elseif x<=3 a=-x;
    else a=-3;
    end

a=0.0;

end

function delta=delta(x)

    if x<=-1 delta=0;
    elseif x<=3 delta=0.5*(x+1);
    else delta = 2;
    end

 delta=0.0;

end



%stiffness for P1
function K=K_lpl(x)
    L=length(x);
    diagon=zeros(L,1);
    diagon_und=zeros(L-1,1);
    diagon(1)=1/(x(2)-x(1));

    diagon(2:1:end-1)=1./(x(2:1:end-1)-x(1:1:end-2))+1./(x(3:1:end)-x(2:1:end-1));
    diagon_und(1:1:end)=-1./(x(2:1:end)-x(1:1:end-1));

    diagon(end-1)=1./(x(end)-x(end-1));

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
        nq=5;
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


function K=construct_block_matrix(x,nu,lambda)
    L=length(x);
    K=zeros(2*L,2*L);
    Mdelta=mass_matrixP1P1_scaled(x, @delta, 3);
    Malpha=mass_matrixP1P1_scaled(x,@alpha,3);
    M=mass_matrixP1(x);
    S=K_lpl(x);
    K(1:L,1:L)=-Malpha-1i*nu*M;
    K(1:L,(L+1):(2*L))=-1i*Mdelta;
    K((L+1):(2*L),1:L)=1i*Mdelta;
    K((L+1):(2*L),(L+1):(2*L))=S-Malpha-1i*nu*M;
    %adding a boundary condition

    K(L+1,L+1)=K(L+1,L+1)-1i*lambda;
end

function b=construct_rhs(x)
    b=zeros(2*length(x),1);
    b(length(x)+1)=exp(1i*x(1));
end






















