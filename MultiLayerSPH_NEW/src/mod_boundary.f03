module boundary_conditions
  use global_var
  implicit none

  ! temporary arrays (temp = BZ->domain, temp2 = domain->BZ)
  real*8 ::           tempV(100), tempV2(100)
  real*8 ::           tempw(100), tempw2(100)
  real*8 ::           temphave(100), temphave2(100)
  real*8 ::           templ(100), templ2(100)
  real*8 ::           tempx(100), tempx2(100)
  real*8 ::           tempa(100), tempa2(100)
  real*8 ::           tempu(100), tempu2(100)
  real*8 ::           tempB(100), tempB2(100)
  real*8 ::           tempinhave(100), tempinhave2(100)
  real*8 ::           tempinl(100), tempinl2(100)
  real*8 ::           temppave(100), temppave2(100)
  integer ::          temploc(100), temploc2(100)
  ! array of empty index
  integer,allocatable ::          idx_p(:)
  integer,allocatable ::          idx_bp(:)

contains
!-------------------- OPEN BOUNDARY --------------------------------------------------------------------
  subroutine open_boundary(p,bp,gr)
    implicit none
    integer ::                          i,k,mm
    type(particles),allocatable ::      bp(:)
    type(particles),allocatable ::      p(:)
    type(mesh),allocatable ::           gr(:)

    do k =1,nlayer
      mm = size(bp(k)%V)
      do i = 1,mm
        if(bp(k)%x(i).lt.domain(1).and.bp(k)%V(i).ne.0.) then
          bp(k)%V(i) = inVsx
          bp(k)%have(i) =hsx/nlayer ! gr(k)%have(idx_init_dom)!
          bp(k)%w(i) =bp(k)%V(i)/bp(k)%have(i)
          bp(k)%u(i) =  inusx(k)*(1d0 + 0.2d0*sin(pi*dble(steptime)*dt))
          !if(steptime*dt >= 0.5d0 .and. steptime*dt <= 0.7d0) bp(k)%u(i) =  (1d0+1d0*cos(2.5*pi*dble(steptime)*dt))*inusx(k) !inusx(k)*(1d0 + 0.2d0*sin(2d0*pi*dble(steptime)*dt))
        
          bp(k)%inhave(i) = bp(k)%have(i)
        else if(bp(k)%x(i).gt.domain(2).and.bp(k)%V(i).ne.0.) then
          ! bp(k)%have(i) = hdx/nlayer!gr(k)%have(idx_end_dom)!
          bp(k)%u(i) = inusx(k)!gr(k)%u(idx_end_dom)!
        end if
      end do
    end do
  end subroutine open_boundary

  subroutine change_particles(p,bp,gr)
    implicit none
    integer ::                          i,k,nn,mm,ii,jj,j_idxp,j_idxbp,ll,pp,sz_p,sz_bp,sz_temp2,sz_temp
    type(particles),allocatable ::      p(:)
    type(particles),allocatable ::      bp(:)
    type(mesh),allocatable ::           gr(:)

    do k = 1,nlayer
      nn = size(p(k)%x)
      mm = size(bp(k)%x)
      ii = 0                    ! index of temporary arrays temp
      ! temporary arrays are set to zero
      tempV = 0
      tempw = 0
      temphave = 0
      templ = 0
      tempx = 0
      tempu = 0
      tempa = 0
      tempB = 0
      tempinhave = 0
      tempinl = 0
      temppave = 0
      temploc = 0

      ! particelle entranti/uscenti nelle buffer zones (in->BZ/BZ->out); per le particelle entranti nella BZ di inflow, devo
      ! necessariamente anche capire quali particelle entrano nel dominio; tali particelle vengono messe nell'array temp
      do i = 1,mm
        if(bp(k)%x(i).gt.gr(k)%x(ninterp) .or. bp(k)%x(i).lt.gr(k)%x(1) .and. bp(k)%V(i).ne.0.) then
          bp(k)%V(i) = 0.
          bp(k)%w(i) = 0.
          bp(k)%have(i) = 0.
          bp(k)%l(i) = 0.
          bp(k)%x(i) = 0.
          bp(k)%u(i) = 0.
          bp(k)%a(i) = 0.
          bp(k)%B(i) = 0.
          bp(k)%inhave(i) = 0.
          bp(k)%inl(i) = 0.
          bp(k)%pave(i) = 0.
          bp(k)%loc(i) = 0.
        else if(bp(k)%x(i).ge.domain(1) .and. bp(k)%x(i).le.domain(2) .and. bp(k)%V(i).ne.0.) then
          ii = ii + 1
          tempV(ii) = bp(k)%V(i)
          tempw(ii) = bp(k)%w(i)
          temphave(ii) = bp(k)%have(i)
          templ(ii) = bp(k)%l(i)
          tempx(ii) = bp(k)%x(i)
          tempu(ii) = bp(k)%u(i)
          tempa(ii) = bp(k)%a(i)
          tempB(ii) = bp(k)%B(i)
          tempinhave(ii) = bp(k)%inhave(i)
          tempinl(ii) = bp(k)%inl(i)
          temppave(ii) = bp(k)%pave(i)
          temploc(ii) = bp(k)%loc(i)

          bp(k)%V(i) = inVsx
          bp(k)%w(i) = 0. ! da definire con bc
          bp(k)%have(i) = 0. ! da definire con bc
          bp(k)%l(i) = gr(k)%l(idx_init_dom)! !gr(k)%l(idx_init_dom)
          bp(k)%x(i) = gr(k)%x(1) + (bp(k)%x(i)-domain(1))
          bp(k)%u(i) = 0. ! da definire con bc
          bp(k)%a(i) = 0. ! calcolate nella subroutine acceleration
          bp(k)%B(i) = 0. ! da definire con gradcorr
          bp(k)%inhave(i) = 0. ! da definire con bc
          bp(k)%inl(i) =  gr(k)%l(1)
          bp(k)%pave(i) = 0. ! da calcolare
          bp(k)%loc(i) = 0 !da calcolare
        end if
      end do

      ! particelle che escono dal dominio e vanno nella BZ di outflow; tali particelle vengono inizialmente messe negli array
      ! temp2, dovranno poi essere spostate negli array bp
      tempV2 = 0.
      tempw2 = 0.
      temphave2 = 0.
      templ2 = 0.
      tempx2 = 0.
      tempu2 = 0.
      tempa2 = 0.
      tempB2 = 0.
      tempinhave2 = 0.
      tempinl2 = 0.
      temploc2 = 0.
      do i = 1,nn
        jj = 0        !index of temporary array temp2
        if(p(k)%x(i).gt.domain(2) .and. p(k)%V(i).ne.0.) then
          jj = jj + 1
          tempV2(jj) = p(k)%V(i)
          tempw2(jj) = p(k)%w(i)
          temphave2(jj) = p(k)%have(i)
          templ2(jj) = p(k)%l(i)
          tempx2(jj) = p(k)%x(i)
          tempu2(jj) = p(k)%u(i)
          tempa2(jj) = p(k)%a(i)
          tempB2(jj) = p(k)%B(i)
          tempinhave2(jj) = p(k)%inhave(i)
          tempinl2(jj) = p(k)%inl(i)
          temploc2(jj) = p(k)%loc(i)

          p(k)%V(i) = 0.
          p(k)%w(i) = 0.
          p(k)%have(i) = 0.
          p(k)%l(i) = 0.
          p(k)%x(i) = 0.
          p(k)%u(i) = 0.
          p(k)%a(i) = 0.
          p(k)%B(i) = 0.
          p(k)%inhave(i) = 0.
          p(k)%inl(i) = 0.
          p(k)%loc(i) = 0
          p(k)%xold(i) = 0.
          p(k)%uold(i) = 0.
          p(k)%aold(i) = 0.
          ! devo mettere a 0 anche pave,pplus,pmin ecc ecc (righe 47-48 mod_global.f03)?
        end if
      end do

      call empty_index(p,bp,k)
      ll = 0
      sz_p = size(idx_p)
      sz_bp = size(idx_bp)
      sz_temp = size(tempV)
      ! trasferisco da temp a p; devo usare idx_p
      do j_idxp = 1,sz_p
        ll = ll + 1
        if(ll.gt.sz_temp) exit

        p(k)%V(idx_p(j_idxp)) = tempV(ll)
        p(k)%w(idx_p(j_idxp)) = tempw(ll)
        p(k)%have(idx_p(j_idxp)) = temphave(ll)
        p(k)%l(idx_p(j_idxp)) = templ(ll)
        p(k)%x(idx_p(j_idxp)) = tempx(ll)
        p(k)%u(idx_p(j_idxp)) = tempu(ll)
        p(k)%B(idx_p(j_idxp)) = tempB(ll)
        p(k)%inhave(idx_p(j_idxp)) = tempinhave(ll)
        p(k)%inl(idx_p(j_idxp)) = tempinl(ll)
        p(k)%pave(idx_p(j_idxp)) = temppave(ll)
        p(k)%loc(idx_p(j_idxp)) = temploc(ll)
        p(k)%a(idx_p(j_idxp)) = tempa(ll)
        p(k)%xold(idx_p(j_idxp)) = tempx(ll)
        p(k)%uold(idx_p(j_idxp)) = tempu(ll)
        p(k)%aold(idx_p(j_idxp)) = tempa(ll)
      end do
      ! trasferisco da p a bp
      pp = 0
      sz_temp2 = size(tempV2)
      do j_idxbp = 1,sz_bp
        pp = pp + 1
        if(pp.gt.sz_temp2) exit
        bp(k)%V(idx_bp(j_idxbp)) = tempV2(pp)
        bp(k)%w(idx_bp(j_idxbp)) = tempw2(pp)
        bp(k)%have(idx_bp(j_idxbp)) = temphave2(pp)
        bp(k)%l(idx_bp(j_idxbp)) = templ2(pp)
        bp(k)%x(idx_bp(j_idxbp)) = tempx2(pp)
        bp(k)%u(idx_bp(j_idxbp)) = tempu2(pp)
        bp(k)%a(idx_bp(j_idxbp)) = tempa2(pp)
        bp(k)%B(idx_bp(j_idxbp)) = tempB2(pp)
        bp(k)%inhave(idx_bp(j_idxbp)) = tempinhave2(pp)
        bp(k)%inl(idx_bp(j_idxbp)) = tempinl2(pp)
        bp(k)%pave(idx_bp(j_idxbp)) = temppave2(pp)
        bp(k)%loc(idx_bp(j_idxbp)) = temploc2(pp)
      end do
    end do
  end subroutine change_particles

  subroutine empty_index(p,bp,k)
    implicit none
    integer ::                          i,k,nn,mm,j_idxp,j_idxbp
    integer,allocatable ::              tempp(:),tempbp(:)
    type(particles),allocatable ::      p(:)
    type(particles),allocatable ::      bp(:)

    if(allocated(idx_p) .and. allocated(idx_bp)) deallocate(idx_p,idx_bp)
    nn = size(p(k)%x)
    mm = size(bp(k)%x)
    allocate(idx_p(nn),idx_bp(mm))
    j_idxp = 0
    j_idxbp = 0
    idx_p = 0
    idx_bp = 0
    ! idx_p
    do i = 1,nn
      if (p(k)%V(i).eq.0.) then
        j_idxp = j_idxp + 1
        idx_p(j_idxp) = i
      end if
    end do
    !idx_bp
    do i = 1,mm
      if (bp(k)%V(i).eq.0.) then
        j_idxbp = j_idxbp + 1
        idx_bp(j_idxbp) = i
      end if
    end do
    allocate(tempp(j_idxp),tempbp(j_idxbp))
    tempp = idx_p(1:j_idxp)
    tempbp =  idx_bp(1:j_idxbp)
    deallocate(idx_p,idx_bp)
    allocate(idx_p(1:j_idxp),idx_bp(1:j_idxbp))
    idx_p(1:j_idxp) = tempp(1:j_idxp)
    idx_bp(1:j_idxbp) = tempbp(1:j_idxbp)
    deallocate(tempp,tempbp)
  end subroutine empty_index

!------------------------------------------------------------------------------------------------------------------

! ---------------------------- CLOSED BOUNDARY --------------------------------------------------------------------
	SUBROUTINE close_boundary(p,bp)
		IMPLICIT NONE
		INTEGER :: i, contatore, contatore2
		INTEGER :: k,nn                                                        ! layer
		DOUBLE PRECISION :: mxv
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

		! the structure "p" is already allocated;
		! in order to allocate "bp" it is necessary to compute the number of particles near the boundary
		do k =1,nlayer
			nn = size(p(k)%x)
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
			DO i = 1,nn
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
			DO i = 1, nn
				IF ((p(k)%x(i)-domain(1)).LE.mxv.AND. p(k)%x(i).NE.domain(1)) THEN
					bp(k)%x(contatore2) = domain(1) - p(k)%x(i)
					bp(k)%u(contatore2) = -p(k)%u(i)
					bp(k)%V(contatore2) = p(k)%V(i)
					bp(k)%w(contatore2) = p(k)%w(i)
					bp(k)%Have(contatore2) = p(k)%Have(i)
					bp(k)%l(contatore2) = p(k)%l(i)
					bp(k)%B(contatore2) = p(k)%B(i)
					bp(k)%Pave(contatore2) = p(k)%Pave(i)
					contatore2 = contatore2 + 1
				ELSE IF ((domain(2)-p(k)%x(i)).LE.mxv .AND. p(k)%x(i).NE.domain(2)) THEN
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
	END SUBROUTINE close_boundary

  SUBROUTINE boundarycupd(p,bp)
		IMPLICIT NONE
		INTEGER :: k, contatore2, i,nn
		DOUBLE PRECISION :: mxv
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
		TYPE(particles), DIMENSION(:), ALLOCATABLE :: bp

		do k = 1,nlayer
			contatore2 = 1
			mxv = 2d0*MAXVAL(p(k)%l)
			nn = size(p(k)%x)

			DO i = 1, nn
				IF ((p(k)%x(i)-domain(1)).LE.mxv.AND. p(k)%x(i).NE.domain(1)) THEN
				  bp(k)%Pave(contatore2) = p(k)%Pave(i)
				  bp(k)%B(contatore2) = p(k)%B(i)
				  contatore2 = contatore2 + 1
				ELSE IF ((domain(2)-p(k)%x(i)).LE.mxv .AND. p(k)%x(i).NE.domain(2)) THEN
				  bp(k)%Pave(contatore2) = p(k)%Pave(i)
				  bp(k)%B(contatore2) = p(k)%B(i)
				  contatore2 = contatore2 + 1
				END IF
			END DO
		end do
  END SUBROUTINE boundarycupd

  SUBROUTINE collision(p)
	IMPLICIT NONE
	INTEGER :: i, k,nn
	DOUBLE PRECISION :: told
	TYPE(particles), DIMENSION(:), ALLOCATABLE :: p
	DO k =1,nlayer
		nn = size(p(k)%x)
		DO i =1,nn
			IF (p(k)%x(i).LT.domain(1)) THEN
				told = ABS(p(k)%x(i))/ABS(p(k)%u(i))
				p(k)%x(i) = told*ABS(p(k)%u(i))
				p(k)%u(i) = -p(k)%u(i)
			ELSE IF (p(k)%x(i).GT.domain(2)) THEN
				told = (p(k)%x(i)-domain(2))/ABS(p(k)%u(i))
				p(k)%x(i) = domain(2)-told*ABS(p(k)%u(i))
				p(k)%u(i) = -p(k)%u(i)
			END IF
		END DO
	END DO
	END SUBROUTINE collision
end module boundary_conditions
