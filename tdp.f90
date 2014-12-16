!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this file contains the program to compute the time dependant pb for Xmode eq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module tdp
implicit none 
contains
  subroutine Kcoeff(x,K1,K2,K1x,K2x,nu,dt,e,me,eps0,B0)
    real(8), intent(in)  :: x,nu,dt,e,me,eps0,B0
    real(8), intent(out) :: K1 , K2 , K1x , K2x 
    
    K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me*eps0)
    K2     =  1 - nu * dt/2 - dt*dt*x*e*e/(4*me*eps0)
    K1x    =  K1 + dt*dt*e*e*B0*B0/(4*K1*me*me)
    K2x    =  K2 -  dt*dt*e*e*B0*B0/(4*K1*me*me)
  end subroutine Kcoeff
  
  subroutine writeall(t,N,Ntime,X,X12,ux,uy,Ex,Ey,H)
    integer,intent (in) :: N, Ntime
    integer, intent(in) :: t
    real(8), dimension(0: N),  intent(in) :: X, ux, Ex, H
    real(8), dimension(0:N-1), intent(in) :: X12, uy, Ey
    character (len=90) :: filename
    integer :: i
    write (filename, '( "password/ux", i7.7, ".data" )' ) t
    open(10, file = filename)
    write (filename, '( "password/Ex", i7.7, ".data" )' ) t
    open(11, file= filename)
    write (filename, '( "password/H", i7.7, ".data" )' ) t
    open(12, file=filename)
    write (filename, '( "password/uy", i7.7, ".data" )' ) t
    open(13, file = filename)
    write (filename, '( "password/Ey", i7.7, ".data" )' ) t
    open(14, file = filename)
    do i = 0, N-1
       write(10,*) X(i), ux(i)
       write(11,*) X(i), Ex(i)
       write(12,*) X(i), H(i)
       write(13,*) X12(i), uy(i)
       write(14,*) X12(i), Ey(i)
    end do
    write(10,*) X(N), ux(N)
    write(11,*) X(N), Ex(N)
    write(12,*) X(N), H(N)
  
  
    
  end subroutine writeall

  
  subroutine tdp_sub(N,Nney,Nhi,Nuxi,Nuy,Neyi,Nex,Nxx,Nx12,dx,dt,Ntime,me,e,eps0,nu,B0,omega,Nei,Ney,uxi,Hi,Exi,Eyi,uyi,X,X12) 
    integer :: N ,Nney, Nhi,Nuxi,Nuy,Neyi,Nex,Nxx,Nx12,Ntime
    real(8) :: dx , dt, me , e , eps0 , nu , B0
    real(8) :: omega
    real(8) :: Nei(0:N),Hi(0:Nhi),Exi(0:Nex),X(0:Nxx)
    real(8) :: uxi(0:Nuxi)
    real(8) :: Ney(0:Nney),Eyi(0:Neyi),uyi(0:Nuy),X12(0:Nx12)
    !f2py optional , depend(Nei)  :: N=len(Nei)
    !f2py optional , depend(Ney)  :: Nney=len(Ney)
    !f2py optional , depend(Hi)   :: Nhi=len(Hi) 
    !f2py optional , depend(uxi)  :: Nuxi=len(uxi)
    !f2py optional , depend(uyi)  :: Nuy=len(uyi)  
    !f2py optional , depend(eyi)  :: Neyi=len(eyi) 
    !f2py optional , depend(exi)  :: Nex=len(exi)  
    !f2py optional , depend(X)    :: Nxx=len(X)
    !f2py optional , depend(X12)  :: Nx12=len(X12)   
    real(8),dimension(0:N) :: ux, tux,H,H12
    real(8) :: Ex(0:N), tEx(0:N)
    real(8) :: Ey(0:N-1),tEy(0:N-1)
    real(8) :: uy(0:N-1),tuy(0:N-1)

    real(8) :: K1, K2, K1x, K2x, Ec, corr,t, En
    complex(16) :: EyF(0:N-1), ExF(0:N-1)  
    

    integer :: i,iter
    print *, "nu", nu
    print *, "dt",dt
    print *, "omega",omega
    t = 0.
    ux(0:N)   = uxi(0:N)
    Ex(0:N)   = Exi(0:N)
    H(0:N)    = Hi(0:N)
    uy(0:N-1) = uyi(0:N-1)
    Ey(0:N-1) = Eyi(0:N-1)
    
    do i=0,N-1
     EyF(i)=cmplx(0.0,0.0);
     ExF(i)=cmplx(0.0,0.0);
    end do

    open(15, file = "password/ET.data")
    open(22, file = "password/ExF.data")
    open(23, file = "password/EyF.data")
    !time loop

    do iter = 0, Ntime
       t = (iter+1.0)*dt
       
       H12(0) = sin(omega*(t-dt/2.0))*(exp(-100*(t-dt/2-1)*(t-dt/2-1)))
       do i = 1, N-1
          H12(i) = H(i) - (dt/dx)*(Ey(i)-Ey(i-1))
       end do
       H12(N) = 0;!H(N)
       
       do i = 0, N-1
          call Kcoeff(NEy(i),K1,K2,K1x,K2x,nu,dt,e,me,eps0,B0)
          tux(i) = (1/K1x) * (K2x*ux(i) +(e/me)*dt*Ex(i) + (((e*e*B0/(me*me))*dt*dt)/(2*K1))*Ey(i)&
               - (e*e*B0/(me*me*eps0))*dt*dt*dt/(4*K1)*((H12(i+1) - H12(i))/dx)   &
               + (e*B0*dt/(2*me))*(K2/K1 + 1)*uy(i))
          
       end do
       call Kcoeff(NEi(N),K1,K2,K1x,K2x,nu,dt,e,me,eps0,B0)
       tux(N) =  (1/K1x) * (K2x*ux(N) +(e/me)*dt*Ex(N)) !&
            !- (e*e*B0/(me*me*eps0))*dt*dt*dt/(4*K1)*(( - H12(N))/dx)  )
       
       do i = 0, N-1
          
          call Kcoeff(NEy(i),K1,K2,K1x,K2x,nu,dt,e,me,eps0,B0)
          tuy(i) = (1/K1) * (K2 * uy(i) &
               +dt * (e/me)* Ey(i) &
               -dt*(dt/2)*(e/(eps0*me)) &
               * (H12(i+1) - H12(i))/dx)&
               -(e*B0/me)*dt*(tux(i)+ux(i))/2
          
       end do
       do i = 0, N-1
          tEx(i) = Ex(i) - dt*(e*NEy(i)/eps0)* (tux(i) + ux(i))/(2)
          tEy(i) = Ey(i) - (dt/eps0) * (H12(i+1) - H12(i))/dx  -(dt/2)*(e*NEy(i)/eps0)*(tuy(i) + uy(i))
       end do
       tEx(N) = Ex(N) - dt*(e*NEi(N)/eps0)* (tux(N) + ux(N))/(2)
       
       ux(0:N) = tux(0:N)
       Ex(0:N) = tEx(0:N)
       uy(0:N-1) = tuy(0:N-1)
       Ey(0:N-1) = tEy(0:N-1)
       H(0:N)  = H12(0:N)
       Ec = 0
       do i = 0, N-1
          Ec = Ec+NEy(i)*(ux(i)*ux(i)+uy(i)*uy(i))
       end do
       corr = 0
       do i = 0, N-2
          corr = corr + (Ey(i+1)-Ey(i))*H12(i+1)
       end do
       En = (dt/dx)*corr + Ec
       do i = 0, N-1
          En = En+ Ey(i)*Ey(i)+Ex(i)*Ex(i) + H12(i)*H12(i)
       end do
       !if (iter>=2e5) then
       call writeall(iter,N,Ntime,X,X12,ux,uy,Ex,Ey,H)            
       !end if
       write(15,*) iter*dt, En
       if (maxval(H)<1E-10) then
          print *, iter,'H nul ?!'
       end if
       
       do i=0, N-1
         EyF(i)=EyF(i)+exp(cmplx(0,1)*omega*t)*Ey(i);
         ExF(i)=ExF(i)+exp(cmplx(0,1)*omega*t)*Ex(i);
       end do
    end do
    call writeall (iter,N,Ntime,X,X12,ux,uy,Ex,Ey,H)
    do i=0,N-1
     write(22,*) ExF(i)
     write(23,*) EyF(i)
    end do
  end subroutine tdp_sub

end module tdp
