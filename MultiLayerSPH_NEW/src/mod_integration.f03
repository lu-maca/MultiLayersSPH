module integration
  use global_var
  implicit none
contains

  subroutine forwardeuler(p,bp)
    implicit none
    integer ::                        k,i, nn

    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)

    !$OMP DO
    do k=1,nlayer
      p(k)%xold = p(k)%x
      p(k)%uold = p(k)%u
      p(k)%u = p(k)%uold + p(k)%a*dt
      p(k)%x = p(k)%xold + p(k)%u*dt

      p(k)%aold = p(k)%a

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

    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)


    DO k=1,nlayer
      p(k)%xold = p(k)%x
      p(k)%uold = p(k)%u
      p(k)%aold = p(k)%a
      
     
      p(k)%u = p(k)%uold + 0.5D0 * (p(k)%a+p(k)%aold)*dt
      p(k)%x = p(k)%xold + p(k)%u*dt + 0.5*p(k)%aold*dt**2
      
      p(k)%aold = p(k)%a

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
          if (p(k)%x(i).le.(domain(1)+1d-5).and.p(k)%uold(i).eq.0.) then
            p(k)%x(i) = domain(1)
            p(k)%u(i) = 0.
          else if (p(k)%x(i).ge.(domain(2)-1d-5).and.p(k)%uold(i).eq.0.) then
            p(k)%x(i) = domain(2)
            p(k)%u(i) = 0.
          end if
    		end do
    	end if
    END DO


  end subroutine leapfrog

end module integration
