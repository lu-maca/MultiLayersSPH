module grid
  use global_var
  use sph_function
  implicit none

contains
  subroutine avepressure(gr)
      implicit none
      integer ::                              i,k
      type(mesh),allocatable ::               gr(:)

      do i = 1,ninterp
        do k =1,nlayer
          gr(k)%pave(i) = 0.5*rho*g*(2.*gr(nlayer)%zcoord(i)-gr(k)%zcoord(i)-gr(k-1)%zcoord(i))*cos(beta)
        end do
      end do
  end subroutine avepressure

  subroutine derivatives(gr)
      implicit none
      integer ::              i,k
      real*8 ::               dx
      real*8,allocatable ::   uave(:)

      type(mesh),allocatable :: gr(:)

      dx = gr(1)%x(2)-gr(1)%x(1)

      do k = 0,nlayer-1
        do i = 1,ninterp
          gr(k+1)%zcoord(i) = gr(k)%zcoord(i)  + gr(k+1)%have(i)
        end do
      end do
      do k = 1,nlayer
        do i = 1,ninterp
          ! CALCOLO DERIVATA
          if(i.eq.1) THEN
            gr(k)%dhdx(i) =(gr(k)%zcoord(i+1)-gr(k)%zcoord(i))/dx
          else if(i.eq.size(gr(k)%zcoord)) then
            gr(k)%dhdx(i) =(gr(k)%zcoord(i)-gr(k)%zcoord(i-1))/dx
          else
            gr(k)%dhdx(i) =0.5D0*(gr(k)%zcoord(i+1)-gr(k)%zcoord(i-1))/dx
          end if
        end do
      end do

      do i = 1,ninterp
        allocate(uave(0:nlayer))
        uave(0) = 0.
        do k = 0,nlayer-1
          uave(k+1) = gr(k+1)%u(i)
        end do
        do k =0,nlayer-1
          if(k.EQ.0) then
         		gr(k)%dudz(i) = 2.*uave(k+1)/gr(k+1)%have(i)
          else if(k.NE.0 .AND. k.NE.nlayer) THEN
            gr(k)%dudz(i) = 2d0*(uave(k+1)-uave(k))/(gr(k+1)%have(i)+gr(k)%have(i))
          end if
        end do
        gr(nlayer)%dudz(i) = 0                                          ! no wind shear
        deallocate(uave)
      end do
  END SUBROUTINE derivatives


! --------------- INITIALIZE GRID -------------------------------
	SUBROUTINE init_grid_CLOSE(gr)
		IMPLICIT NONE
		DOUBLE PRECISION :: dx
		INTEGER :: i, k
		TYPE(mesh), DIMENSION(:), ALLOCATABLE :: gr

		dx = (domain(2)-domain(1))/(ninterp-1)
		DO k = 0,nlayer
	    IF (steptime .NE. 0) THEN
	        DEALLOCATE(gr(k)%x,gr(k)%u,gr(k)%have,gr(k)%dudz,gr(k)%dhdx,gr(k)%l,gr(k)%pint,gr(k)%pave)
	    END IF
	    ALLOCATE(gr(k)%x(ninterp))
	    ALLOCATE(gr(k)%u(ninterp))
	    ALLOCATE(gr(k)%have(ninterp))
	    ALLOCATE(gr(k)%dudz(ninterp))
	    ALLOCATE(gr(k)%dhdx(ninterp))
	    ALLOCATE(gr(k)%l(ninterp))
	    ALLOCATE(gr(k)%pint(ninterp))
	    ALLOCATE(gr(k)%pave(ninterp))
	    ALLOCATE(gr(k)%zcoord(ninterp))
		END DO

		DO k = 0,nlayer
	    DO i = 1,ninterp
    	  gr(k)%x(i) = domain(1) + (i-1)*dx
	    END DO
	    gr(k)%u = 0
	    gr(k)%have = 0
	    gr(k)%dudz = 0
	    gr(k)%dhdx = 0
	    gr(k)%l = 0
	    gr(k)%pint = 0
	    gr(k)%pave = 0
	    gr(k)%zcoord = 0
	END DO
 	END SUBROUTINE init_grid_CLOSE


  subroutine init_grid_OPEN(gr)
    implicit none
    integer ::                        ni,i,k,np
    real*8 ::                         dxpart,dxgr,Lbz
    type(mesh),allocatable ::         gr(:)

    np = numpar
    dxpart = (domain(2)-domain(1))/(np+1)
    Lbz = 4.*dxpart

    ni = ninterp

    do k = 0,nlayer
      allocate(gr(k)%x(ni),gr(k)%u(ni),gr(k)%have(ni),gr(k)%dudz(ni))
      allocate(gr(k)%dhdx(ni),gr(k)%l(ni),gr(k)%pint(ni),gr(k)%pave(ni))
      allocate(gr(k)%zcoord(ni))
    end do

    do k = 0,nlayer
      dxgr = (domain(2)-domain(1)+2.*Lbz)/(ni-1)
      do i = 1,ni
          gr(k)%x(i) = domain(1)-Lbz + (i-1)*dxgr
      end do
      gr(k)%u = 0
      gr(k)%have = 0
      gr(k)%dudz = 0
      gr(k)%dhdx = 0
      gr(k)%l = 0
      gr(k)%pint = 0
      gr(k)%pave = 0
      gr(k)%zcoord = 0
    end do

    do i = 1,ninterp
      if(gr(1)%x(i).ge.domain(1)) then
        idx_end_dom = i-1
        exit
      endif
    end do
    do i = 1,ninterp
      if(gr(1)%x(ninterp-i).le.domain(2)) then
        idx_init_dom = i-1
        exit
      endif
    end do
  end subroutine init_grid_OPEN

! --------------------------------------------------------------

  subroutine interp(p,bp,gr)
    implicit none
    integer ::                        sxi, i, k, nn,mm
    real*8 ::                         griddx, gridsx, parpos
    type(particles),allocatable ::    p(:)
    type(particles),allocatable ::    bp(:)
    type(mesh),allocatable ::         gr(:)


    do k = 1,nlayer
      nn = size(p(k)%x)
      do i = 1,nn
        if(p(k)%V(i).ne.0.) then
          sxi = p(k)%loc(i)
          gridsx = gr(k)%x(sxi)
          griddx = gr(k)%x(sxi+1)
          parpos = p(k)%x(i)

          p(k)%Pave(i) = gr(k)%pave(sxi)+(gr(k)%pave(sxi+1)-gr(k)%pave(sxi))/(griddx-gridsx)*(parpos-gridsx)
          p(k)%pplus(i) = gr(k)%pint(sxi)+(gr(k)%pint(sxi+1)-gr(k)%pint(sxi))/(griddx-gridsx)*(parpos-gridsx)
          p(k)%dhdxplus(i) = gr(k)%dhdx(sxi)+(gr(k)%dhdx(sxi+1)-gr(k)%dhdx(sxi))/(griddx-gridsx)*(parpos-gridsx)
          p(k)%dudzplus(i) = gr(k)%dudz(sxi)+(gr(k)%dudz(sxi+1)-gr(k)%dudz(sxi))/(griddx-gridsx)*(parpos-gridsx)

          p(k)%pmin(i) = gr(k-1)%pint(sxi)+(gr(k-1)%pint(sxi+1)-gr(k-1)%pint(sxi))/(griddx-gridsx)*(parpos-gridsx)
          p(k)%dhdxmin(i) = gr(k-1)%dhdx(sxi)+(gr(k-1)%dhdx(sxi+1)- gr(k-1)%dhdx(sxi))/(griddx-gridsx)*(parpos-gridsx)
          p(k)%dudzmin(i) = gr(k-1)%dudz(sxi)+(gr(k-1)%dudz(sxi+1)-gr(k-1)%dudz(sxi))/(griddx-gridsx)*(parpos-gridsx)
        end if
      end do

      if (boundary.eqv. .false.) then
		    nn = size(bp(k)%x)
		    do i = 1, nn
		      if(bp(k)%V(i).ne.0.) then
		        sxi = bp(k)%loc(i)
		        gridsx = gr(k)%x(sxi)
		        griddx = gr(k)%x(sxi+1)
		        parpos = bp(k)%x(i)
		        bp(k)%Pave(i) = gr(k)%pave(sxi)+(gr(k)%pave(sxi+1)-gr(k)%pave(sxi))/(griddx-gridsx)*(parpos-gridsx)
        	end if
     		end do
      end if
    end do
  end subroutine interp

  subroutine locpar(p,bp,gr)
    implicit none
    integer ::                              i,l,k,nn
    real*8 ::                               valore,valstep
    type(mesh),allocatable ::               gr(:)
    type(particles),allocatable ::          p(:)
    type(particles),allocatable ::          bp(:)

    valstep = gr(1)%x(2)-gr(1)%x(1)
    !internal particles
    do k = 1, nlayer
      nn = size(p(k)%x)
      do i = 1,nn
        valore = p(k)%x(i)
        do l = 1,ninterp
          if (boundary.eqv. .true. .and. abs(valore-gr(1)%x(l)).le. valstep .and. valore.eq.domain(2)) then
            p(k)%loc(i) = l-1
          else if(abs(valore-gr(1)%x(l)).le. valstep .and. (valore-gr(1)%x(l)).ge.0) then
            p(k)%loc(i) = l
            exit
          else if(ABS(valore-gr(1)%x(l)).le. valstep .and. (valore-gr(1)%x(l)).lt.0) then
            p(k)%loc(i) = l-1
            exit
          end if
        end do
      end do
    end do

		if (boundary.eqv. .false.) then
		  ! external particles
		  do k = 1, nlayer
		    nn = size(bp(k)%x)
		    do i = 1,nn
		      valore = bp(k)%x(i)
		      do l = 1,ninterp
		        if(abs(valore-gr(1)%x(l)).le. valstep .and. (valore-gr(1)%x(l)).ge.0) then
		          bp(k)%loc(i) = l
		          if (bp(k)%loc(i).eq.ninterp) bp(k)%loc(i)=ninterp-1
		          exit
		        else if(ABS(valore-gr(1)%x(l)).le. valstep .and. (valore-gr(1)%x(l)).LT.0) then
		          bp(k)%loc(i) = l-1
		          exit
		        end if
		      end do
		    end do
		  end do
		end if
  end subroutine locpar

  subroutine pressure(gr)
    implicit none
    integer ::                    i, k
    type(mesh),allocatable ::     gr(:)

    ! CALCOLO PRESSIONE
    gr(nlayer)%pint = patm
    do k =nlayer-1,0,-1
        do i = 1, ninterp
            gr(k)%pint(i) = gr(k+1)%pint(i) + rho*g*gr(k+1)%have(i)*cos(beta)
        end do
    end do
  end subroutine pressure

! ------------- COMPUTE QUANTITIES AT GRID NODES ----------------

  subroutine quantitiesatgridnodes_OPEN(p,bp,gr)
    implicit none
    integer ::              i,k,m,loc,nn,mm
    real*8 ::               u,tl,R,Rbc,sommad
    real*8,allocatable ::   Rvec(:)
    type(mesh),allocatable ::               gr(:)
    type(particles),allocatable ::          p(:)
    type(particles),allocatable ::          bp(:)

    do k = 1,nlayer
      nn = size(p(k)%x)
      mm = size(bp(k)%x)
      do i = 1,ninterp
        if (gr(k)%x(i).ge.domain(1).and.gr(k)%x(i).le.domain(2)) then
          if (allocated(Rvec)) deallocate(Rvec)
          allocate(Rvec(nn))
          Rvec = huge(1.)
          do m = 1,nn
            if(p(k)%V(m).ne.0.) Rvec(m) = abs(gr(k)%x(i) - p(k)%x(m))
          end do
          loc = minloc(Rvec,dim = 1)
          gr(k)%l(i) = p(k)%l(loc)
        else if(gr(k)%x(i).lt.domain(1)) then
          if (allocated(Rvec)) deallocate(Rvec)
          allocate(Rvec(size(bp(k)%x)))
          Rvec = huge(1.)
          do m = 1,mm
            if(bp(k)%V(m).ne.0.) Rvec(m) = abs(gr(k)%x(i) - bp(k)%x(m))
          end do
          loc = minloc(Rvec,dim = 1)
          gr(k)%l(i) = bp(k)%l(loc)
        else if(gr(k)%x(i).gt.domain(2)) then
          if (allocated(Rvec)) deallocate(Rvec)
          allocate(Rvec(mm))
          Rvec = huge(1.)
          do m = 1,mm
            if(bp(k)%V(m).ne.0.) Rvec(m) = abs(gr(k)%x(i) - bp(k)%x(m))
          end do
          loc = minloc(Rvec,dim = 1)
          gr(k)%l(i) = bp(k)%l(loc)
        end if
      end do
    end do

    do k = 1,nlayer
      nn = size(p(k)%x)
      mm = size(bp(k)%x)
      do i =1,ninterp
        if (gr(k)%x(i).ge.domain(1).and.gr(k)%x(i).le.domain(2)) then
          u = 0
          tl = 0
          R = 0
          Rbc = 0
          do m = 1,nn
              R = gr(k)%x(i) - p(k)%x(m)
              if (ABS(R) .LE. 2.*gr(k)%l(i).and. p(k)%V(m).ne.0.) THEN
                  u = u + p(k)%w(m)*p(k)%u(m)*Ww(R,gr(k)%l(i))
                  tl = tl + p(k)%w(m)*p(k)%Have(m)*Ww(R,gr(k)%l(i))
              end if
          end do
          do m = 1, mm
              Rbc = gr(k)%x(i) - bp(k)%x(m)
              if (ABS(Rbc).LE.2.*gr(k)%l(i).and. bp(k)%V(m).ne.0.) THEN
                  u = u + bp(k)%w(m)*bp(k)%u(m)*Ww(Rbc,gr(k)%l(i))
                  tl = tl + bp(k)%w(m)*bp(k)%Have(m)*Ww(Rbc,gr(k)%l(i))
              end if
          end do
          gr(k)%u(i) = u
          gr(k)%have(i) = tl
        else if(gr(k)%x(i).lt.domain(1)) then
          sommad = 0
          u = 0
          tl = 0
          R = 0
          Rbc = 0
          do m = 1,nn
              R = gr(k)%x(i) - p(k)%x(m)
              if (abs(R) .le. 2.*gr(k)%l(i) .and. p(k)%V(m).ne.0.) then
                  u = u + p(k)%w(m)*p(k)%u(m)*Ww(R,gr(k)%l(i))
                  tl = tl + p(k)%w(m)*p(k)%Have(m)*Ww(R,gr(k)%l(i))
                  sommad = sommad + p(k)%w(m)*Ww(R,gr(k)%l(i))
              end if
          end do
          do m = 1, mm
              Rbc = gr(k)%x(i) - bp(k)%x(m)
              if (abs(Rbc).le.2.*gr(k)%l(i).and. bp(k)%V(m).ne.0.) then
                  u = u + bp(k)%w(m)*bp(k)%u(m)*Ww(Rbc,gr(k)%l(i))
                  tl = tl + bp(k)%w(m)*bp(k)%Have(m)*Ww(Rbc,gr(k)%l(i))
                  sommad = sommad + bp(k)%w(m)*Ww(Rbc,gr(k)%l(i))
              end if
          end do
          gr(k)%u(i) = u/sommad
          gr(k)%have(i) = tl/sommad
        else if (gr(k)%x(i).gt.domain(2)) then
          sommad = 0
          u = 0
          tl = 0
          R = 0
          Rbc = 0
          do m = 1,nn
              R = gr(k)%x(i) - p(k)%x(m)
              if (abs(R) .le. 2.*gr(k)%l(i) .and. p(k)%V(m).ne.0.) then
                  u = u + p(k)%w(m)*p(k)%u(m)*Ww(R,gr(k)%l(i))
                  tl = tl + p(k)%w(m)*p(k)%Have(m)*Ww(R,gr(k)%l(i))
                  sommad = sommad + p(k)%w(m)*Ww(R,gr(k)%l(i))
              end if
          end do
          do m = 1, mm
              Rbc = gr(k)%x(i) - bp(k)%x(m)
              if (abs(Rbc).le.2.*gr(k)%l(i).and. bp(k)%V(m).ne.0.) then
                  u = u + bp(k)%w(m)*bp(k)%u(m)*Ww(Rbc,gr(k)%l(i))
                  tl = tl + bp(k)%w(m)*bp(k)%Have(m)*Ww(Rbc,gr(k)%l(i))
                  sommad = sommad + bp(k)%w(m)*Ww(Rbc,gr(k)%l(i))
              end if
          end do
          gr(k)%u(i) = u/sommad
          gr(k)%have(i) = tl/sommad
          !gr(k)%u(i) = gr(k)%u(idx_end_dom)
          !gr(k)%have(i) = gr(k)%have(idx_end_dom)
        end if
      end do
    end do
  end subroutine quantitiesatgridnodes_OPEN


	SUBROUTINE quantitiesatgridnodes_CLOSE(p,bp,gr)
		IMPLICIT NONE
		INTEGER :: i,m,k, loc
		DOUBLE PRECISION :: R, Rbc, uu,tl
		DOUBLE PRECISION, DIMENSION(Np) :: Rvec
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp
		TYPE(mesh), DIMENSION(:), ALLOCATABLE :: gr

		do k =1,nlayer
			Rvec = 0
			DO i = 1,ninterp
				DO m = 1,Np
					Rvec(m) = ABS(gr(k)%x(i) - p(k)%x(m))
				END DO
				loc = MINLOC(Rvec, DIM = 1)
				gr(k)%l(i) = p(k)%l(loc)
			END DO
			DO i = 1,ninterp
				! typical sph approx
				uu = 0
				tl = 0
				R = 0
				Rbc = 0
				DO m = 1,Np
					R = gr(k)%x(i) - p(k)%x(m)
					IF (ABS(R) .LE. 2d0*gr(k)%l(i)) THEN
						uu = uu + p(k)%w(m)*p(k)%u(m)*Ww(R,gr(k)%l(i))
						tl = tl + p(k)%w(m)*p(k)%Have(m)*Ww(R,gr(k)%l(i))
					END IF
				END DO
				DO m = 1, size(bp(k)%x)
					Rbc = gr(k)%x(i) - bp(k)%x(m)
					IF (ABS(Rbc).LE.2d0*gr(k)%l(i)) THEN
						uu = uu + bp(k)%w(m)*bp(k)%u(m)*Ww(Rbc,gr(k)%l(i))
						tl = tl + bp(k)%w(m)*bp(k)%Have(m)*Ww(Rbc,gr(k)%l(i))
					END IF
				END DO
				gr(k)%u(i) = uu
				gr(k)%have(i) = tl
			END DO
		end do
	END SUBROUTINE quantitiesatgridnodes_CLOSE

end module grid
