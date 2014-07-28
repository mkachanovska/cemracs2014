function [u,K,x]=solve_m(dx,lambda, nu)
%mesh generation
L=2;
H=10;

x=-L:dx:H;


%phi=linspace(-pi/2.0,pi/2.0, 100);
%x=(H+L)*sin(phi)/2.0+(H-L)./2.0;

K=full_stiffness(x,lambda,nu);
r=rhs(x);
disp( "formed the matrix and the RHS");
u=K\r;
u=[u; 0];
end


function a=alpha(x)

if x<=-1 a=1;
elseif x<=3 a=-x;
else a=-3;
end

end

function delta=delta(x)

if x<=-1 delta=0;
elseif x<=3 delta=0.5*(x+1);
else delta = 2;
end

end

function f=perturb(x,nu)

if x<=-1 f=-(1+1i*nu);
elseif x<=3 f=1/4*(x+1)^2*(x-1i*nu)./(x^2+nu^2)-x-1i*nu;
else f=4*(-3-1i*nu)./(9+nu^2)-3-1i*nu;

end



function K=K_lpl(x)
L=length(x);
diagon=zeros(L-1,1);
diagon_und=zeros(L-2,1);
diagon(1)=1/(x(2)-x(1));
for k=2:1:L-1
diagon(k)=1/(x(k)-x(k-1))+1/(x(k+1)-x(k));
diagon_und(k-1)=-1/(x(k)-x(k-1));
end

K=diag(diagon)+diag(diagon_und,-1)+diag(diagon_und,1);
end

%returns the part related to delta, alpha





function K=K_MassScaled(x, nu)
%unautomated procedure, for these particular values of delta and alpha



end

function K=K_general(x,nu)
    diagon=zeros(length(x)-1,1);
    subdiagon=zeros(length(x)-2,1);
    updiagon=zeros(length(x)-1,1);
end



function K=K_const(x,nu)

end




function K=full_stiffness(x, lambda,nu)
%K=K_NonConstQuad(x,nu)+K_lpl(x);
K=K_lpl(x);
K(1,1)=K(1,1)+1i*lambda;



end




	
function r=rhs(x)
r=zeros(length(x)-1,1);
r(1)=-2;
end






