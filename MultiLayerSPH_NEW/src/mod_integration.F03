module integration
  use global_var
  implicit none
contains

  subroutine CFL_routine(p,gr)
    implicit none
    integer :: i, k, nn(nlayer), nnn
    real*8 :: dt_cfl, dt_mu, diff_stab(nlayer)
    real*8, allocatable :: den(:,:)
    type(particles), allocatable :: p(:)
    type(mesh), allocatable :: gr(:)

    do k = 1,nlayer
      nn(k) = size(p(k)%x)
    end do
    nnn = maxval(nn)
    dt_prev= dt
    ! cfl condition
    if (allocated(den)) deallocate(den)
    allocate(den(nlayer,nnn))
    do k = 1,nlayer
      do i = 1,nnn
        if (p(k)%V(i).ne.0d0) then
          den(k,i) = p(k)%l(i)/(sqrt(g*p(k)%have(i)) + p(k)%u(i))
        else
          den(k,i) = huge(1.)
        end if
      end do
    end do
    dt_cfl = minval(cfl_fixed*den)

    ! diffusion condition
    do k=1,nlayer
      if (nu.eq.0.) then
        diff_stab(k) = huge(1.)
      else
        diff_stab(k) = (minval(gr(k)%have,mask = gr(k)%have.gt.0. ))**2/(2*nu)
      end if
    end do
    dt_mu = minval(diff_stab)

    if (dt_cfl<= dt_mu) then
      dt = dt_cfl
    else
      dt = dt_mu
    end if

    if (time + dt>tfin) then
        dt = abs( -tfin+time-dt+dt_prev)
    endif
    if (dt .le. 0d0) then
      print*,dt,dt_mu,dt_cfl
      stop "dt less than zero"
    end if
  end subroutine CFL_routine


  subroutine forwardeuler(p,bp)
    implicit none
    integer ::                        k,i, nn

    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)

    !$OMP DO
    do k=1,nlayer
      p(k)%aold = p(k)%a
      p(k)%uold = p(k)%u
      p(k)%xold = p(k)%x

      p(k)%u = p(k)%uold + p(k)%a*dt
      p(k)%x = p(k)%x + p(k)%u*dt
      

      if (boundary.eqv..false.) then
        bp(k)%x = bp(k)%x + bp(k)%u*dt + bp(k)%a*dt**2
      else
        nn = size(p(k)%x)
        do i =1,nn
          if (p(k)%xold(i).eq.domain(1)) then
            p(k)%x(i) = domain(1)
            p(k)%u(i) = 0.
          else if (p(k)%xold(i).eq.domain(2)) then
            p(k)%x(i) = domain(2)
            p(k)%u(i) = 0.
          end if
    		end do
      end if
    end do
    !$OMP END DO
  end subroutine forwardeuler

  subroutine leapfrog(p,bp)
    implicit none
    integer ::                        k,i,nn,mm
    real*8, allocatable :: xx(:)

    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)


    DO k=1,nlayer
     
      p(k)%uold = p(k)%u
      p(k)%aold = p(k)%a
!      if (boundary.eqv. .true.) then
!          allocate(xx(size(p(k)%x)))
!          xx = p(k)%x
!          p(k)%x = xx + (xx - p(k)%xold)*dt/dt_prev + 0.5d0*p(k)%a*(dt+dt_prev)*dt
!           p(k)%xold = xx
!       deallocate(xx)
!       else
          p(k)%u = p(k)%uold + p(k)%aold*dt
          p(k)%x = p(k)%x + p(k)%u*dt
!       endif
      if (boundary.eqv..false.) then
        bp(k)%x = bp(k)%x + bp(k)%u*dt + bp(k)%a*dt**2
        nn = size(p(k)%x)
        do i =1,nn
          if ( abs(p(k)%x(i)).ge.huge(1.)) then
          print*, i,k,p(k)%x(i), p(k)%u(i),p(k)%aold(i)
          stop 'Error integration line 60'
          end if
        end do
      else
        nn = size(p(k)%x)
				do i =1,nn
          if (p(k)%x(i).le.(domain(1)+1d-3).and.p(k)%uold(i).eq.0.) then
            p(k)%x(i) = domain(1)
            p(k)%u(i) = 0.
          else if (p(k)%x(i).ge.(domain(2)-1d-3).and.p(k)%uold(i).eq.0.) then
            p(k)%x(i) = domain(2)
            p(k)%u(i) = 0.
          end if
    		end do
    	end if
    END DO


  end subroutine leapfrog

end module integration
