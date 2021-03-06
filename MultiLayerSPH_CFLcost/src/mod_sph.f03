module sph_function
  use global_var
  implicit none
contains

  subroutine acceleration(p,bp)
    implicit none
    real*8 ::                         S, A, R
    real*8 ::                         brackdx, bracksx, n1,n2,n3,n4, visc
    integer ::                        i,j,k,nn,mm
    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)

    do k = 1,nlayer
    	nn=size(p(k)%x)
    	mm = size(bp(k)%x)
      do i = 1,nn
        if(p(k)%V(i).ne.0d0) then
          S = g*sin(beta)*p(k)%Have(i)+&
              &1d0/rho*(p(k)%pplus(i)*p(k)%dhdxplus(i)- p(k)%pmin(i)*p(k)%dhdxmin(i))+&
              &nu*(p(k)%dudzplus(i)-p(k)%dudzmin(i))
          A = 0
          do j =1,nn
            if(p(k)%V(j).ne.0.) then
              R = p(k)%x(i)-p(k)%x(j)
              if(abs(R).le.2d0*p(k)%l(i)) then
                  n4 = (p(k)%u(i)-p(k)%u(j))*(p(k)%x(i)-p(k)%x(j))
                  if(n4.GE.0d0) then
                      visc = 0d0
                  else
                      n1 = 0.5d0*(p(k)%l(i)+p(k)%l(j))
                      n2 = 0.5d0*(p(k)%Have(i) + p(k)%Have(j))
                      n3 = 0.5d0*(sqrt(g*p(k)%Have(i)) + sqrt(g*p(k)%Have(j)))
                      visc = -art_visc_coeff*n1*n2*n3*n4/(p(k)%x(i)-p(k)%x(j))**2
                  end if
                  bracksx = 1d0/rho*(p(k)%Have(i))*p(k)%Pave(i) + 0.5d0*visc
                  brackdx = 1d0/rho*(p(k)%Have(j))*p(k)%Pave(j) + 0.5d0*visc
                  A = A + p(k)%w(j)*(bracksx*p(k)%B(i)*dWw(R,p(k)%l(i)) + brackdx*p(k)%B(j)*dWw(R,p(k)%l(j)))
              end if
            end if
          end do
          do j = 1,mm
            if(bp(k)%V(j).ne.0.) then
              R = (p(k)%x(i)-bp(k)%x(j))
              if(abs(R).le.2d0*p(k)%l(i)) then
                  n4 = (p(k)%u(i)-bp(k)%u(j))*(p(k)%x(i)-bp(k)%x(j))
                  if (n4.ge.0d0) then
                      visc = 0d0
                  else
                      n1 = 0.5d0*(p(k)%l(i) + bp(k)%l(j))
                      n2 = 0.5d0*(p(k)%Have(i) + bp(k)%Have(j))
                      n3 = 0.5d0*(sqrt(g*p(k)%Have(i)) + sqrt(g*bp(k)%Have(j)))
                      visc = -art_visc_coeff*n1*n2*n3*n4/(p(k)%x(i)-bp(k)%x(j))**2
                  end if
                  bracksx = 1d0/rho*p(k)%Have(i)*p(k)%Pave(i) + 0.5d0*visc
                  brackdx = 1d0/rho*bp(k)%Have(j)*bp(k)%Pave(j) + 0.5d0*visc
                  A = A + bp(k)%w(j)*(bracksx*p(k)%B(i)*dWw(R,p(k)%l(i)) + brackdx*bp(k)%B(j)*dWw(R,bp(k)%l(j)))
              end if
            end if
          end do
          p(k)%a(i) =1d0/p(k)%have(i)*(S-A)

          if (abs(p(k)%a(i)) .ge. huge(1.)) then
            print*,'acceleration inf',p(k)%have(i),S,A,p(k)%x(i),p(k)%xold(i)
            stop
          elseif (isnan(p(k)%a(i))) then
            print*, 'acceleration nan at', p(k)%x(i)
          endif
        end if
      end do

      if (boundary .eqv. .false.) then
        do i = 1,mm
          if (bp(k)%V(i).ne.0.) bp(k)%a(i) = 0d0
        end do
      endif
    end do
  end subroutine acceleration

  subroutine secant(p,bp)
    IMPLICIT NONE
      INTEGER :: i,k, j,mm,nn
      DOUBLE PRECISION :: old, oold, fold, foold, now, lold, loold,infty,sommaoldd, sommaooldd
      DOUBLE PRECISION :: sommaold,sommaoold,eps,toll,r
      TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
      TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp


      infty = huge(1.)
      DO k =1,nlayer
        mm = size(bp(k)%x)
		nn = size(p(k)%x)
         DO i = 1,nn
            if(p(k)%V(i).ne.0.) then
           ! secant method
             eps = 1D-15
             toll = 10d0
             sommaold = 0d0
             sommaoldd = 0d0
             sommaoold = 0d0
             sommaooldd = 0d0
             lold = 0d0
             oold = p(k)%have(i)
             loold = p(k)%l(i)

             old = p(k)%have(i)*(1d0 - 0.5d0)

             lold = p(k)%inl(i)*p(k)%inhave(i)/old

             foold = 0d0
             fold = 0d0
             DO j = 1,nn
               r = p(k)%x(i)-p(k)%x(j)
                IF (abs(r).le.2d0*p(k)%l(i).and. p(k)%V(j).ne.0.) then
                  sommaooldd = sommaooldd + p(k)%w(j)*Ww(r,loold)
                  sommaold = sommaold + p(k)%V(j)*Ww(r,lold)
                end if
              END DO
              DO j = 1,mm
                r = p(k)%x(i)-bp(k)%x(j)
                IF (abs(r).le.2d0*p(k)%l(i).and. bp(k)%V(j).ne.0.) then
                  sommaoold = sommaoold + bp(k)%V(j)*Ww(r,loold)
                  sommaold = sommaold + bp(k)%V(j)*Ww(r,lold)
                end if
              END DO
              fold = old - sommaold
              foold = oold - sommaoold
              now = old - fold*(old-oold)/(fold-foold)
              oold = old
              old = now
              loold = lold
              lold = p(k)%inl(i)*p(k)%inhave(i)/old

              DO WHILE (abs(toll).ge. eps)
                sommaold = 0d0
                sommaoold = 0d0
                DO j = 1,nn
                  r = p(k)%x(i)-p(k)%x(j)
                  IF (abs(r).le.2d0*p(k)%l(i).and. p(k)%V(j).ne.0.) then
                    sommaoold = sommaoold + p(k)%V(j)*Ww(r,loold)
                    sommaold = sommaold + p(k)%V(j)*Ww(r,lold)
                  end if
                END DO
                DO j = 1,mm
                  r = p(k)%x(i)-bp(k)%x(j)
                  IF (abs(r).le.2d0*p(k)%l(i).and. bp(k)%V(j).ne.0.) then
                    sommaoold = sommaoold + bp(k)%V(j)*Ww(r,loold)
                    sommaold = sommaold + bp(k)%V(j)*Ww(r,lold)
                  end if
                END DO
                fold = old - sommaold
                foold = oold - sommaoold

                now = old - fold*(old-oold)/(fold-foold)
                oold = old
                old = now
                loold = lold
                lold = p(k)%inl(i)*p(k)%inhave(i)/old

                toll = (old-oold)/old
                p(k)%have(i) = old
                p(k)%l(i) = lold

                if (p(k)%l(i).ge.huge(1.)) then
                  print*,'posizione', p(k)%x(i)
                  stop 'l inf'
                elseif (p(k)%l(i).le. 0.) then
                  print*, p(k)%l(i),old
                  stop 'l negative'
                elseif (isnan(p(k)%l(i))) then
                  print*, "l nan", i,k
                  stop "l nan"
                endif
                p(k)%w(i) = p(K)%V(i)/p(k)%have(i)


              END DO
          end if
        END DO
      END DO

  end subroutine secant

  subroutine update_w(p,bp)
    integer :: i, k, j,nn,mm
    real*8 :: somma,R
    TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
    TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

    do k=1,nlayer
      nn = size(p(k)%x)
			mm = size(bp(k)%x)
      do i = 1,nn
        if (p(k)%V(i).ne.0.) then
          somma = 0d0
          do j=1,nn
            R = p(k)%x(i) - p(k)%x(j)
            IF (ABS(R).LE.2D0*p(k)%l(i)) THEN
              somma = somma + p(k)%w(j)*(p(k)%w(j) - p(k)%w(i))*dWw(R,p(k)%l(i))
            END IF
          end do
          do j=1,mm
            R = p(k)%x(i) - bp(k)%x(j)
            IF (ABS(R).LE.2D0*p(k)%l(i).and. bp(k)%V(j).ne. 0.) THEN
              somma = somma + bp(k)%w(j)*(bp(k)%w(j) - p(k)%w(i))*dWw(R,p(k)%l(i))
            END IF
          end do

          p(k)%w(i) = p(k)%w(i)*exp(somma*dt)
          p(k)%have(i) = p(k)%V(i)/p(k)%w(i)
        endif
      end do
    end do
  end subroutine update_w

! ------------------- GRADIENT CORRECTION ---------------------------------

	SUBROUTINE gradient_correction_CLOSE(p,bp)
		IMPLICIT NONE
		INTEGER :: i,j,k,nn,mm
		DOUBLE PRECISION :: bb, R, infinity
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

		! calcolo B per le particelle interne
		infinity = huge(1d0)
		DO k = 1,nlayer
			nn = size(p(k)%x)
			mm = SIZE(bp(k)%x)
			DO i = 1, nn
				bb = 0d0
				DO j = 1,nn
					R = p(k)%x(i) - p(k)%x(j)
					IF (ABS(R).LE.2D0*p(k)%l(i)) THEN
						bb = bb + p(k)%w(j)*(p(k)%x(j)-p(k)%x(i))*dWw(R,p(k)%l(i))
					END IF
				END DO
				DO j = 1, mm
					R = p(k)%x(i) - bp(k)%x(j)
					IF (ABS(R).LE.2D0*p(k)%l(i)) THEN
						bb = bb + bp(k)%w(j)*(bp(k)%x(j)-p(k)%x(i))*dWw(R,p(k)%l(i))
					END IF
				END DO
				p(k)%B(i) = 1D0/bb
				if (p(k)%B(i).GT.infinity) then
					print*,"grad correction ","part",i,"layer",k, "time",steptime,"l", p(k)%l(i)
          print*, p(k)%l
					stop
				endif
			END DO
		END DO
	END SUBROUTINE gradient_correction_CLOSE

  subroutine gradient_correction_OPEN(p,bp)
    implicit none
    integer ::        i,j,k,nn,mm
    real*8 ::         infty,r,bden,specpos
    type(particles),allocatable ::        p(:)
    type(particles),allocatable ::        bp(:)

    infty = huge(1d0)
    ! internal particles
    do k = 1,nlayer
      nn = size(p(k)%x)
      mm =size(bp(k)%x)
      do i = 1,nn
        bden = 0
        do j = 1,nn
          if(p(k)%V(i).ne.0. .and.p(k)%V(j).ne.0.) then
            r = p(k)%x(i) - p(k)%x(j)
            if(abs(r).le.2.*p(k)%l(i)) then
              bden = bden + p(k)%w(j)*(p(k)%x(j)-p(k)%x(i))*dWw(R,p(k)%l(i))
            end if
          end if
        end do
        do j = 1,mm
          if(p(k)%V(i).ne.0. .and. bp(k)%V(j).ne.0.) then
            r = p(k)%x(i) - bp(k)%x(j)
            if(abs(R).le.2.*p(k)%l(i)) then
              bden = bden + bp(k)%w(j)*(bp(k)%x(j)-p(k)%x(i))*dWw(R,p(k)%l(i))
            end if
          end if
        end do
      !  if (bden.eq.0.) continue
        if(p(k)%V(i).ne.0.) then
          p(k)%B(i) = 1./bden
        end if
        if (p(k)%B(i).gt.infty) then
          print*,'p',i,k,bden,steptime
          print*, "particle position (real - boundary)"
          print*,p(k)%x
          print*, bp(k)%x
          print*, "particle volume (real - boundary)"
          print*,p(k)%V
          print*, bp(k)%V
          print*, "smoothing length (real - boundary)"
          print*, p(k)%l
          print*, bp(k)%l
          stop "subroutine gradient_correction_OPEN"
        endif
      end do
    end do
	  !external particles

	  do k = 1,nlayer
	    do i =1,mm
	      bden = 0
	      if(bp(k)%V(i).ne.0.) then
	        if(bp(k)%x(i).lt.domain(1)) then
	          specpos = domain(1)-bp(k)%x(i)
	        else if(bp(k)%x(i).gt.domain(2)) then
	          specpos = 2.*domain(2)-bp(k)%x(i)
	        end if
	      else
	        cycle
	      end if

	      do j = 1,nn
	        if(bp(k)%V(i).ne.0. .and. p(k)%V(j).ne.0.) then
	          r = specpos - p(k)%x(j)
	          if(abs(r).le.2.*bp(k)%l(i)) then
	            bden = bden + p(k)%w(j)*(p(k)%x(j)-specpos)*dWw(R,bp(k)%l(i))
	          !  print*,bden,i,k
	          end if
	        end if
	      end do
	      do j = 1,mm
	        if(bp(k)%V(i).ne.0. .and. bp(k)%V(j).ne.0.) then
	          r = specpos - bp(k)%x(j)
	          if(abs(r).le.2.*bp(k)%l(i)) then
	            bden = bden + bp(k)%w(j)*(bp(k)%x(j)-specpos)*dWw(R,bp(k)%l(i))
	          !  print*,bden,i,k
	          end if
	        end if
	      end do
        !if (bden.eq.0.) continue
	      if(p(k)%V(i).ne.0.) then
	        bp(k)%B(i) = 1./bden
	      end if

	      if (bp(k)%B(i).gt.infty) then
	        print*,'bp',i,k,bden,steptime
	        print*,bp(k)%x
	        stop
	      endif
	    end do
	  end do
  end subroutine gradient_correction_OPEN

	!------------------------------------------------------------------------------------

	! ----------- initialize particles (position, velocity and properties) ---------------
  SUBROUTINE set_particles_CLOSE(p,bp)
    IMPLICIT NONE
    DOUBLE PRECISION :: hlayer, dx,htot, dam_pos,xmid,amp_pert,sigma
    INTEGER :: i, k, contatore = 1
    TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
    TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

		Np = numpar
    DO k = 1, nlayer
      ALLOCATE(p(k)%x(Np),p(k)%u(Np),p(k)%xold(Np),p(k)%uold(Np),p(k)%aold(Np),p(k)%a(Np))
      ALLOCATE(p(k)%V(Np),p(k)%w(Np),p(k)%Have(Np),p(k)%l(Np),p(k)%B(Np),p(k)%loc(Np))
      ALLOCATE(p(k)%Pave(Np),p(k)%pplus(Np),p(k)%pmin(Np),p(k)%dudzplus(Np),p(k)%dudzmin(Np))
      ALLOCATE(p(k)%dhdxplus(Np),p(k)%dhdxmin(Np),p(k)%inhave(np),p(k)%inl(np))
      ALLOCATE(bp(k)%x(contatore),bp(k)%u(contatore),bp(k)%V(contatore),bp(k)%w(contatore))
      ALLOCATE(bp(k)%Have(contatore),bp(k)%l(contatore),bp(k)%B(contatore),bp(k)%Pave(Np))
    END DO

    DO k = 1,nlayer
      p(k)%x = 0
      p(k)%u = 0                                                         ! potrebbe essere sensato imporre nusselt
      p(k)%a= 0
      p(k)%B = 1
      p(k)%uold = 0
      p(k)%xold = 0
      p(k)%aold = 0

      select case (testcase)
      case (1) ! dam break with dry bed
        htot = 1d0
        hlayer = htot/nlayer
        dam_pos = 5d0
        dx = dam_pos/(Np-1)
        p(k)%V = hlayer*dx
        p(k)%l = 2.5d0*dx ! coefficient must be > 2.2, otherwise the secant method does not converge
        p(k)%w = dx
      case (2) ! dam break with wet bed
        htot = 1 ! left state with h_L > h_R
        hlayer = htot/nlayer
        dx = (domain(2)-domain(1))/(Np-1)
        dam_pos = 5
        p(k)%V = hlayer*dx
        p(k)%V(int(dam_pos/dx)+1:np) = 0.5*hlayer*dx ! the 0.5 adjusts the height of the right state
        p(k)%l =  2.5*dx
        p(k)%w = dx
        p(k)%w(int(dam_pos/dx)+1:np) = 0.5*dx
      case (3)
        stop "Wrong boundary conditions! You cannot choose a Nusselt flow with wall boundaries. Change input file settings."
      case (4)
        htot = 1 ! undisturbed depth
        sigma = 0.5 ! std dev
        amp_pert = 0.5 ! disturbance amplitude
        xmid = (domain(2)-domain(1))/2.
        dx = (domain(2)-domain(1))/(Np-1)
        hlayer = htot/nlayer
        p(k)%l =  2*dx
        p(k)%w = dx
      case default
        print*, "Choose an existing test case"
      end select
    END DO

    DO k =1,nlayer
      DO i = 1,Np
        p(k)%x(i) = domain(1) + (i-1)*dx
        if (testcase.eq.4) p(k)%V(i) = hlayer*dx*(1 + amp_pert* exp(-((p(k)%x(i) - xmid)/sigma)**2) )

        p(k)%Have(i) = p(k)%V(i)/p(k)%w(i)
        p(k)%inhave(i) = p(k)%Have(i)
        p(k)%inl(i) = p(k)%l(i)
      END DO
    END DO
  END SUBROUTINE set_particles_CLOSE

  subroutine set_particles_OPEN(p,bp)
    implicit none
    integer ::                              k,i,np
    real*8 ::                               hlayer,dx,htot
    type(particles),allocatable ::          p(:)
    type(particles),allocatable ::          bp(:)

    do k = 1,nlayer
      p(k)%np = numpar
    end do

    ! allocate arrays
    do k = 1,nlayer
      np = p(k)%np
      ! fluid particles
      allocate(p(k)%x(np+20),p(k)%u(np+20),p(k)%a(np+20))
      allocate(p(k)%V(np+20),p(k)%w(np+20),p(k)%have(np+20),p(k)%l(np+20))
      allocate(p(k)%B(np+20),p(k)%inhave(np+20),p(k)%inl(np+20))
      allocate(p(k)%xold(np+20),p(k)%uold(np+20),p(k)%aold(np+20))
      allocate(p(k)%pave(np+20),p(k)%pplus(np+20),p(k)%pmin(np+20))
      allocate(p(k)%dudzplus(np+20),p(k)%dudzmin(np+20),p(k)%loc(np+20))
      allocate(p(k)%dhdxplus(np+20),p(k)%dhdxmin(np+20))

      !open boundary particles
      allocate(bp(k)%x(np),bp(k)%u(np),bp(k)%a(np))
      allocate(bp(k)%V(np),bp(k)%w(np),bp(k)%have(np),bp(k)%l(np))
      allocate(bp(k)%B(np),bp(k)%inhave(np),bp(k)%inl(np))
      allocate(bp(k)%pave(np),bp(k)%loc(np))
      allocate(bp(k)%xold(np),bp(k)%uold(np),bp(k)%aold(np)) !non dovrebbero servire
    end do

    do k = 1,nlayer
        dx = (domain(2)-domain(1))/(np+1)
        p(k)%x = 0d0
        p(k)%u = 0d0
        p(k)%a= 0d0
        p(k)%B = 1d0
        p(k)%uold = 0d0
        p(k)%xold = 0d0
        p(k)%aold = 0d0
    end do

    hsx = 0
    hdx = 0
    allocate(inusx(nlayer),inudx(nlayer))
    do k =1,nlayer
        do i = 1,np
          select case (testcase)
          case (1)
            stop "Wrong boundary conditions!"
          case (2)
            stop "Wrong boundary conditions!"
          case (3) ! nusselt
            if (beta.eq. 0 .or. nu.eq.0.) stop "Nusselt flow requires non null viscosity and angle of inclination:&
              & check if nu and/or beta are set equal to zero in the input file."
            htot = 0.05d0
            hlayer = htot/nlayer
            p(k)%V(i) = hlayer*dx

            p(k)%w(i) = dx
            p(k)%x(i) =  dx + domain(1) + (i-1)*dx
            p(k)%l(i) =  2d0*dx
            p(k)%u(i) =  .5d0*g/nu*(2*htot-(dble(k)-0.5)*hlayer)*(dble(k)-0.5)*hlayer*sin(beta)
            p(k)%have(i) = p(k)%V(i)/p(k)%w(i)
            p(k)%inhave(i) = p(k)%Have(i)
            p(k)%inl(i) = p(k)%l(i)
          case (4)
            stop "Wrong boundary conditions!"
          end select
        end do
        hsx = hsx + p(k)%inhave(1)
        hdx = hdx + p(k)%inhave(p(k)%Np)

        inusx(k) = p(k)%u(1)
        inudx(k) = p(k)%u(p(1)%np)
    end do
    inVsx = p(1)%V(1)
    inVdx = p(1)%V(p(1)%np)
  end subroutine set_particles_OPEN

  !-----------------------------------------------------------------------------------

  subroutine set_bparticles_OPEN(p,bp,gr)
    implicit none
    integer ::          i,k,nn
    real*8 ::           Lbz
    type(mesh),allocatable ::               gr(:)
    type(particles),allocatable ::          p(:)
    type(particles),allocatable ::          bp(:)

    do k = 1,nlayer
      Lbz = domain(1)-gr(k)%x(1)
      nn = size(p(k)%x)
      ! particella al contorno sx
      bp(k)%x(nn/2) = domain(1)
      bp(k)%u(nn/2) = p(k)%u(1)
      bp(k)%V(nn/2) = p(k)%V(1)
      bp(k)%w(nn/2) = p(k)%w(1)
      bp(k)%have(nn/2) = p(k)%have(1)
      bp(k)%l(nn/2) = p(k)%l(1)
      bp(k)%B(nn/2) = p(k)%B(1)
      bp(k)%inhave(nn/2) = p(k)%inhave(1)
      bp(k)%inl(nn/2) = p(k)%inl(1)
      ! particella al contorno dx
      bp(k)%x(nn/2+1) = domain(2)
      bp(k)%u(nn/2+1) = p(k)%u(1)
      bp(k)%V(nn/2+1) = p(k)%V(1)
      bp(k)%w(nn/2+1) = p(k)%w(1)
      bp(k)%have(nn/2+1) = p(k)%have(1)
      bp(k)%l(nn/2+1) = p(k)%l(1)
      bp(k)%B(nn/2+1) = p(k)%B(1)
      bp(k)%inhave(nn/2+1) = p(k)%inhave(1)
      bp(k)%inl(nn/2+1) = p(k)%inl(1)
      do i = 1,nn
        if ((p(k)%x(i)-domain(1)).lt.Lbz .and. p(k)%x(i).ne.domain(1)) then
          bp(k)%x(i) = domain(1)-p(k)%x(i)
          bp(k)%u(i) = p(k)%u(i)
          bp(k)%V(i) = p(k)%V(i)
          bp(k)%w(i) = p(k)%w(i)
          bp(k)%have(i) = p(k)%have(i)
          bp(k)%l(i) = p(k)%l(i)
          bp(k)%B(i) = p(k)%B(i)
          bp(k)%inhave(i) = p(k)%inhave(i)
          bp(k)%inl(i) = p(k)%inl(i)
        else if (abs(p(k)%x(i)-domain(2)).lt.Lbz .and. p(k)%x(i).ne.domain(2)) then
          bp(k)%x(i) = 2.*domain(2)-p(k)%x(i)
          bp(k)%u(i) = p(k)%u(i)
          bp(k)%V(i) = p(k)%V(i)
          bp(k)%w(i) = p(k)%w(i)
          bp(k)%have(i) = p(k)%have(i)
          bp(k)%l(i) = p(k)%l(i)
          bp(k)%B(i) = p(k)%B(i)
          bp(k)%inhave(i) = p(k)%inhave(i)
          bp(k)%inl(i) = p(k)%inl(i)
        end if
      end do
    end do
  end subroutine set_bparticles_OPEN

	SUBROUTINE set_bparticles_CLOSE(p,bp)
		IMPLICIT NONE
		INTEGER :: i, contatore, contatore2
		INTEGER :: k                                                        ! layer
		DOUBLE PRECISION :: mxv
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

		! the structure "p" is already allocated;
		! in order to allocate "bp" it is necessary to compute the number of particles near the boundary
		do k=1,nlayer
			DEALLOCATE(bp(k)%x)
			DEALLOCATE(bp(k)%u)
			DEALLOCATE(bp(k)%V)
			DEALLOCATE(bp(k)%w)
			DEALLOCATE(bp(k)%Have)
			DEALLOCATE(bp(k)%l)
			DEALLOCATE(bp(k)%B)
			DEALLOCATE(bp(k)%Pave)
			contatore = 0
			mxv = 2d0*MAXVAL(p(k)%l)
			DO i = 1,size(p(k)%x)
				IF ((p(k)%x(i)-domain(1)).LE.mxv .AND. p(k)%x(i).NE.domain(1)) THEN
				  contatore = contatore + 1
				ELSEIF ((domain(2)-p(k)%x(i)).LE.mxv .AND. p(k)%x(i).NE.domain(2)) THEN
				  contatore = contatore + 1
				END IF
			END DO

			ALLOCATE(bp(k)%x(contatore))
			ALLOCATE(bp(k)%u(contatore))
			ALLOCATE(bp(k)%V(contatore))
			ALLOCATE(bp(k)%w(contatore))
			ALLOCATE(bp(k)%Have(contatore))
			ALLOCATE(bp(k)%l(contatore))
			ALLOCATE(bp(k)%B(contatore))
			ALLOCATE(bp(k)%Pave(contatore))

			contatore2 = 1
			DO i = 1,size(p(k)%x)
				IF ((p(k)%x(i)-domain(1)).LE.mxv.AND. p(k)%x(i).NE.domain(1)) THEN !
					bp(k)%x(contatore2) = domain(1) - p(k)%x(i)
					bp(k)%u(contatore2) = -p(k)%u(i)
					bp(k)%V(contatore2) = p(k)%V(i)
					bp(k)%w(contatore2) = p(k)%w(i)
					bp(k)%Have(contatore2) = p(k)%Have(i)
					bp(k)%l(contatore2) = p(k)%l(i)
					bp(k)%B(contatore2) = p(k)%B(i)
					bp(k)%Pave(contatore2) = p(k)%Pave(i)
					contatore2 = contatore2 + 1
				ELSE IF ((domain(2)-p(k)%x(i)).LE.mxv .AND. p(k)%x(i).NE.domain(2)) THEN !
					bp(k)%x(contatore2) = 2D0*domain(2) - p(k)%x(i)
					bp(k)%u(contatore2) = -p(k)%u(i)
					bp(k)%V(contatore2) = p(k)%V(i)
					bp(k)%w(contatore2) = p(k)%w(i)
					bp(k)%Have(contatore2) = p(k)%Have(i)
					bp(k)%l(contatore2) = p(k)%l(i)
					bp(k)%B(contatore2) = p(k)%B(i)
					bp(k)%Pave(contatore2) = p(k)%Pave(i)
					contatore2 = contatore2 + 1
				END IF
			END DO
		end do
	END SUBROUTINE set_bparticles_CLOSE

!--------------------------------------------------------------------

  subroutine output(gr,p,bp)
    implicit none
    character(len=10) ::    file_id
    character(len=50) ::    file_name
    integer ::              i,k,nn,mm

    type(mesh), allocatable :: gr(:)
    type(particles), allocatable :: p(:)
    type(particles), allocatable :: bp(:)


    ! scrittura file
    write(file_id, '(i0)') steptime
    file_name = 'file' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(file_name), status='new')

    do i = 1,ninterp
      write(1,*) gr(1)%x(i), gr(nlayer)%zcoord(i),gr(1)%zcoord(i)
    end do
    close(1)

    if (mod(steptime,1*res_freq).eq.0) THEN
      write(file_id, '(i0)') steptime
      file_name = 'vel' // trim(adjustl(file_id)) // '.dat'
      open(2, file=trim(file_name),STATUS='new')
      file_name = 'height'// trim(adjustl(file_id)) // '.dat'
      open(3,file=trim(file_name),STATUS='new')
      file_name = 'zcoord'// TRIM(adjustl(file_id)) // '.dat'
      OPEN(4,file=TRIM(file_name),STATUS='new')
!      file_name = 'upart'// TRIM(adjustl(file_id)) // '.dat'
!      OPEN(5,file=TRIM(file_name),STATUS='new')
!      file_name = 'smotlength'// TRIM(adjustl(file_id)) // '.dat'
!      OPEN(50,file=TRIM(file_name),STATUS='new')
!      file_name = 'xpart'// TRIM(adjustl(file_id)) // '.dat'
!      OPEN(30,file=TRIM(file_name),STATUS='new')
!      file_name = 'apart'// TRIM(adjustl(file_id)) // '.dat'
!      OPEN(90,file=TRIM(file_name),STATUS='new')

      do k = 1,nlayer
        write(2,*) gr(k)%u(:)
      end do
      do k =1,nlayer
        write(3,*) gr(k-1)%zcoord(:) + 0.5D0*gr(k)%have(:)
      end do
      DO k =1,nlayer
          WRITE(4,*) gr(k)%zcoord(:)
      END DO

!    k = 1                              ! strato in output
!      nn = size(p(k)%x)
!      mm = size(bp(k)%x)
!    DO i =1,nn
!        WRITE(5,*) p(k)%u(i)
!    END DO
!    DO i = 1,mm
!        WRITE(5,*) bp(k)%u(i)
!    END DO
!    DO i =1,nn
!        WRITE(50,*) p(k)%l(i)
!    END DO
!    DO i = 1,mm
!        WRITE(50,*) bp(k)%l(i)
!    END DO
!    ! WRITE(30,*) "particles"
!     DO i =1,nn
!         WRITE(30,*) p(k)%x(i), p(k)%have(i), p(k)%V(i), p(k)%B(i)
!     END DO
!     ! WRITE(30,*) "boundary particles"
!    DO i = 1,mm
!       WRITE(30,*) bp(k)%x(i), bp(k)%have(i), bp(k)%V(i), bp(k)%B(i)
!    END DO

      close(2)
      close(3)
      close(4)
!      close(5)
!      close(50)
!      close(30)
    end if
  end subroutine output

  ! subroutine stability(p,gr)
  !   implicit none
  !   real*8,allocatable :: cfl(:),den(:),diff_stab(:)
  !   integer :: k, i,j,nn
  !   type(particles), allocatable :: p(:)
  !   type(mesh), allocatable :: gr(:)
  !
  !   allocate(cfl(nlayer),diff_stab(nlayer))
  !   do k = 1,nlayer
  !     nn = size(p(k)%x)
  !     if (allocated(den)) deallocate(den)
  !     allocate(den(nn))
  !     do i = 1,nn
  !       den(i) = p(k)%l(i)/(sqrt(g*p(k)%have(i)) + p(k)%u(i))
  !     end do
  !     cfl(k) = dt/minval(den)
  !     !if (minval(den) == 0) cfl(k) = 0d0
  !
  !     if (nu.eq.0.) then
  !       diff_stab(k) = huge(1.)
  !     else
  !       diff_stab(k) = (minval(gr(k)%have,mask = gr(k)%have.gt.0. ))**2/nu
  !     end if
  !   end do
  !
  !   if (mod(steptime,1*res_freq).eq.0) then
  !     if (steptime .eq.0 ) then
  !       open(200, file = 'CFL.dat',STATUS='new')
  !       write(200,*) 'CFL<1   ', 'diffusion stability <0.5    ', 'timestep    ', 'physical time'
  !       write(200,*) maxval(CFL), dt * 1./minval(diff_stab), steptime, steptime*dt
  !       close(200)
  !     else
  !       open(200, file = 'CFL.dat',STATUS='old', position='append')
  !       write(200,*) maxval(CFL), dt * 1./minval(diff_stab), steptime, steptime*dt
  !       close(200)
  !     end if
  !   end if
  !
  !   if (maxval(cfl).ge.1 .or. dt * 1./minval(diff_stab).ge.0.5 ) then
  !     if (testcase.ne. 1) then
  !       print*, 'Stability conditions not respected: ', 'CFL=',maxval(cfl), 'diffusion=',  dt * 1./minval(diff_stab)
  !       print*, 'WARNING: CFL > 1'
  !     end if
  !   endif
  !
  ! end subroutine

! ---------------- SPH kernel -----------------------------------
  function Ww(r,h)
    real*8 ::        Ww, r, h
    if(abs(r)>=0 .and. abs(r)<=h) then
        Ww = 2/(3*h)*(1-1.5*(abs(r)/h)**2 * (1-abs(r)/(2*h)))
    else if(ABS(r)>h .and. ABS(r)<=2*h) then
        Ww = 1/(6*h)*(2-abs(r)/h)**3
    else if(ABS(r)>2*h) then
        Ww = 0
    end if
  end function Ww

  function dWw(r,h)
    real*8 ::        dWw, r, h
    if(abs(r)>=0 .and. abs(r)<=h) then
      dWw = (r*(-2*h + 1.5*abs(r)))/h**4
    else if(abs(r)>h .and. abs(r)<=2*h) then
      dWw = - r*(abs(r)-2*h)**2/(2*h**4*abs(r))
    else if(abs(r)>2*h) then
      dWw = 0
    end if
  end function dWw

end module sph_function
