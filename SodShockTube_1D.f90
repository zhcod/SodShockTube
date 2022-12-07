! This code is to solve the Sod Shock tube problem (1D).
! Using Lax, Lax-Wendroff, MacC, FVS, Roe method.
! Write by zhliu on 2022.11


module sodshock 
    implicit none 
    integer,parameter           ::  Nx = 1000
    real*8,parameter            ::  xl = 1.D0
    real*8,parameter            ::  end_time = 0.2D0
    real*8,parameter            ::  dx = xl/Nx
    real*8,parameter            ::  CFL=0.8D0
    real*8,parameter            ::  Gamma=1.4D0
    contains

    subroutine SetInitalValue(u)
        implicit none 
        real*8,dimension(3,0:Nx+1)    ::  u
        u(1,0:Nx/2)=1.D0 
        u(1,Nx/2+1:Nx+1)=0.125D0
        u(2,0:Nx+1)=0.D0 
        u(3,0:Nx/2)=1.D0/(Gamma-1.D0)
        u(3,Nx/2+1:Nx+1)=0.1D0/(Gamma-1.D0)
    end subroutine

    subroutine SetBoundaryCondition(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u
        u(:,0)=u(:,1)
        u(:,Nx+1)=u(:,Nx)
    end subroutine

    subroutine SetBoundary(u)
        implicit none
        real*8,dimension(0:Nx+1)    ::  u
        u(0)=u(1)
        u(Nx+1)=u(Nx)
    end subroutine

    subroutine ComputeTimeStep(u,dt)
        implicit none 
        real*8,dimension(3,0:Nx+1)    ::  u
        real*8,dimension(0:Nx+1)      ::  rho,v,p,a
        real*8     ::  dt
        call ComputeOriginalVariable(u,rho,v,p)
        a = sqrt(Gamma*p/rho)
        dt = CFL*dx/maxval(abs(v)+a)
    end subroutine

    subroutine ComputeOriginalVariable(u,rho,v,p)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u
        real*8,dimension(0:Nx+1)      ::  rho,v,p
        integer  :: j
        do j=0,Nx+1 
            rho(j)  =   u(1,j)
            v(j)    =   u(2,j)/u(1,j)
            p(j)    =   (Gamma-1)*(u(3,j)-0.5D0*u(2,j)**2/u(1,j))
        end do 
    end subroutine

    subroutine ComputeFlux(u,f)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,f
        real*8,dimension(0:Nx+1)      ::  rho,v,p
        integer                       ::  j
        call ComputeOriginalVariable(u,rho,v,p)
        do j=0,Nx+1
            f(1,j)=u(2,j)
            f(2,j)=0.5D0*(3.D0-Gamma)*u(2,j)**2/u(1,j)+(Gamma-1.D0)*u(3,j)
            f(3,j)=Gamma*u(2,j)*u(3,j)/u(1,j)-0.5D0*(Gamma-1.D0)*u(2,j)**3/u(1,j)**2
        end do  
    end subroutine

    subroutine Write2File(u,FileIndex,FileName)
        implicit none
        integer       ::   FileIndex
        character*7  ::   FileName  
        integer       ::   j
        real*8,dimension(3,0:Nx+1)    ::  u
        real*8,dimension(0:Nx+1)      ::  rho,v,p
        call ComputeOriginalVariable(u,rho,v,p)
        open(unit=FileIndex, file=FileName)
        write(FileIndex,*) 'rho',',','v',',','p'
        do j=1,Nx 
            write(FileIndex,*) rho(j),',',v(j),',',p(j)
        end do
        close(FileIndex)
    end subroutine

    subroutine Solvewith_LAX(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,f,uu
        integer                     ::  j,k
        real*8                     :: t=0.D0
        real*8  :: dt
        do while(t.le.end_time)
            call ComputeTimeStep(u,dt)
            call ComputeFlux(u,f)
            do j=1,Nx
                do k=1,3
                    uu(k,j)=0.5D0*(u(k,j+1)+u(k,j-1))-0.5D0*(dt/dx)*(f(k,j+1)-f(k,j-1))
                end do 
            end do
            call SetBoundaryCondition(uu)
            u=uu
            t=t+dt 
        end do 
    end subroutine

    subroutine Solvewith_LW(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,f,uh,fh,uu
        integer                     ::  j,k,i
        real*8  :: t=0.D0
        real*8  :: dt
        do while(t.le.end_time)
            call ComputeTimeStep(u,dt)
            call ComputeFlux(u,f)
            do j=1,Nx
                do k=1,3
                    uh(k,j)=0.5D0*(u(k,j+1)+u(k,j))-0.5D0*(dt/dx)*(f(k,j+1)-f(k,j))
                end do 
            end do
            call SetBoundaryCondition(uh)
            call ComputeFlux(uh,fh)
            do j=1,Nx
                do k=1,3
                    uu(k,j)=u(k,j)-(dt/dx)*(fh(k,j)-fh(k,j-1))
                end do 
            end do
            call SetBoundaryCondition(uu)
            u=uu
            t=t+dt
        end do 
    end subroutine

    subroutine Solvewith_MacC(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,f,uh,fh,uu
        integer                     ::  j,k,i
        real*8  :: t=0.D0
        real*8  :: dt 
        do while(t.le.end_time)
            call ComputeTimeStep(u,dt)
            call ComputeFlux(u,f)
            do j=1,Nx
                do k=1,3
                    uh(k,j)=u(k,j)-(dt/dx)*(f(k,j+1)-f(k,j))
                end do 
            end do
            call SetBoundaryCondition(uh)
            call ComputeFlux(uh,fh)
            do j=1,Nx
                do k=1,3
                    uu(k,j)=0.5D0*(u(k,j)+uh(k,j))-0.5D0*(dt/dx)*(fh(k,j)-fh(k,j-1))
                end do 
            end do
            call SetBoundaryCondition(uu)
            u=uu
            t=t+dt
        end do 
    end subroutine

    subroutine Solvewith_FVS(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,uu,f,fp,fm
        real*8,dimension(0:Nx+1)      ::  rho,v,p,a
        integer                     ::  j,k,i
        real*8  :: t=0.D0
        real*8  :: dt 
        do while(t.le.end_time)
            call ComputeTimeStep(u,dt)
            call ComputeOriginalVariable(u,rho,v,p)
            a = sqrt(Gamma*p/rho)
            call ComputeFlux(u,f)

            fp(1,:)=rho/2.D0/Gamma*((2.D0*Gamma-1.D0)*v+a)
            fp(2,:)=rho/2.D0/Gamma*(2.D0*(Gamma-1.D0)*v**2+(v+a)**2)
            fp(3,:)=rho/2.D0/Gamma*((Gamma-1.D0)*v**3+(v+a)**3/2.D0+(3-Gamma)*(v+a)*a**2/(Gamma-1)/2.D0)
            fm(1,:)=rho/2.D0/Gamma*(v-a)
            fm(2,:)=rho/2.D0/Gamma*((v-a)**2)
            fm(3,:)=rho/2.D0/Gamma*((v-a)**3/2.D0+(3-Gamma)*(v-a)*a**2/(Gamma-1)/2.D0)

            do j=1,Nx
                if(v(j).gt.a(j)) then
                write(*,*) 'FVS_ERROR!'
                stop 
                else
                end if 
                do k=1,3
                    uu(k,j)=u(k,j)-(dt/dx)*(fp(k,j)-fp(k,j-1))-(dt/dx)*(fm(k,j+1)-fm(k,j))
                end do 
            end do
            call SetBoundaryCondition(uu)
            u=uu
            t=t+dt
        end do 
    end subroutine


    subroutine Solvewith_Roe(u)
        implicit none
        real*8,dimension(3,0:Nx+1)    ::  u,uu,dh,malpha,ml,mlp,mlm
        real*8,dimension(0:Nx+1)      ::  rho,v,p,a,h,e
        real*8,dimension(0:Nx+1)      ::  mrho,mv,mp,ma,mh,mD
        real*8,dimension(0:Nx+1)      ::  drho,dv,dp
        real*8,dimension(3,3,0:Nx+1)      ::  mR
        integer                     ::  j,k,i
        real*8  :: t=0.D0
        real*8  :: dt 
        do while(t.le.end_time)
            call ComputeTimeStep(u,dt)
            call ComputeOriginalVariable(u,rho,v,p)
            a = sqrt(Gamma*p/rho)
            e = p/(Gamma-1.D0)/rho
            h = e+v**2/2.D0+p/rho

            !前缀m表示Roe平均
            do j=1,Nx
            mD(j)=sqrt(rho(j+1)/rho(j))
            mh(j)=(h(j)+mD(j)*h(j+1))/(1+mD(j))
            mrho(j)=rho(j)*(0.5D0*(1.D0+mD(j)))**2
            mv(j)=(v(j)+mD(j)*v(j+1))/(1+mD(j))
            ma(j)=sqrt((Gamma-1)*(mh(j)-mv(j)**2/2.D0))
            drho(j)=rho(j+1)-rho(j)
            dp(j)=p(j+1)-p(j)
            dv(j)=v(j+1)-v(j)
            end do 
            call SetBoundary(mrho)
            call SetBoundary(mv)
            call SetBoundary(mp)
            call SetBoundary(drho)
            call SetBoundary(dv)
            call SetBoundary(dp)

            malpha(1,:)=drho-dp/ma**2
            malpha(2,:)=0.5D0/ma**2*(dp-mrho*ma*dv)
            malpha(3,:)=0.5D0/ma**2*(dp+mrho*ma*dv)

            ml(1,:)=mv
            ml(2,:)=mv-ma
            ml(3,:)=mv+ma
            mlp=0.5D0*(ml+abs(ml))
            mlm=0.5D0*(ml-abs(ml))

            mR(1,1,:) = 1.D0
            mR(2,1,:) = mv
            mR(3,1,:) = 0.5D0*mv**2

            mR(1,2,:) = 1.D0
            mR(2,2,:) = mv - ma
            mR(3,2,:) = mh - mv*ma

            mR(1,3,:) = 1.D0
            mR(2,3,:) = mv + ma
            mR(3,3,:) = mh + mv*a

            do j=1,Nx
                do k=1,3
                    dh=0.D0
                    do i=1,3
                    dh(k,j)=dh(k,j)+mlm(i,j)*malpha(i,j)*mR(k,i,j)+mlp(i,j-1)*malpha(i,j-1)*mR(k,i,j-1)
                    end do      
                    uu(k,j)=u(k,j)-(dt/dx)*dh(k,j)
                end do 
            end do
            call SetBoundaryCondition(uu)
            u=uu
            t=t+dt
        end do 
    end subroutine

end module 

    
program main 
    use sodshock
    implicit none
    real*8,dimension(3,0:Nx+1)    ::  u

    call SetInitalValue(u)
    call  Solvewith_LAX(u)
    call     Write2File(u,11,'Lax.csv')

    call SetInitalValue(u)
    call   Solvewith_LW(u)
    call     Write2File(u,12,'L_W.csv')

    call SetInitalValue(u)
    call Solvewith_MacC(u)
    call     Write2File(u,13,'Mac.csv')

    call SetInitalValue(u)
    call  Solvewith_FVS(u)
    call     Write2File(u,14,'FVS.csv')

    call SetInitalValue(u)
    call  Solvewith_Roe(u)
    call     Write2File(u,15,'Roe.csv')


end program main





