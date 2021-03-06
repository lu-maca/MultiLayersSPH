program solver
! -------------------------------
! -------------------------------
! Multilayer shallow water solver
! by Luca Macavero
! -------------------------------
! -------------------------------
  use global_var
  use sph_function
  use grid
  use boundary_conditions
  use integration

  implicit none
  real ::                               start,finish
  integer ::                            totstep,k


  type(particles),allocatable ::        p(:)
  type(particles),allocatable ::        bp(:)
  type(mesh),allocatable ::             gr(:)

  call cpu_time(start)

  steptime = 0
  call read_input()
  allocate(gr(0:nlayer),p(nlayer),bp(nlayer))

	select case (boundary)
		case (.true.) !closed
			do k = 1,nlayer
				p(k)%np = Np
			end do
			call set_particles_CLOSE(p,bp)
			call init_grid_CLOSE(gr)
			call set_bparticles_CLOSE(p,bp)

		case (.false.) !open
			call set_particles_OPEN(p,bp)
			call init_grid_OPEN(gr)
			call set_bparticles_OPEN(p,bp,gr)
  end select

  open(unit=10, file = "test.dat")
  write(10,*) "smoothing length iniziale",p(1)%l(1)
  write(10,*) "dt ", dt
  write(10,*) "numero particelle iniziale", p(1)%np
  write(10,*) "nodi", ninterp
  if (SLupd) then
  	write(10,*) "smoothing length update YES"
  else
  	write(10,*) "smoothing length update NO"
  end if
  close(10)
    !call output(gr,p,bp)

  select case (boundary)
	case (.true.)
	  call quantitiesatgridnodes_CLOSE(p,bp,gr)
	  call gradient_correction_CLOSE(p,bp)
	case (.false.)
	  call quantitiesatgridnodes_OPEN(p,bp,gr)
	  call gradient_correction_OPEN(p,bp)
  end select

  call locpar(p,bp,gr)
  call derivatives(gr)
  call pressure(gr)
  call avepressure(gr)
  call output(gr,p,bp)
  call interp(p,bp,gr)

  if (boundary.eqv..true.) then ! closed
  	call boundarycupd(p,bp)!close_boundary(p,bp)
  end if

  call turbulent_viscosity(p,gr)
  call stability(p,gr)
  
  call acceleration(p,bp)
  call forwardeuler(p,bp)

	select case (boundary)
		case (.true.) !closed
			call collision(p)
			call close_boundary(p,bp)
		case (.false.) !open
			call change_particles(p,bp,gr)
		  call open_boundary(p,bp,gr)
	end select

	if (SLupd .eqv. .true.) then
		select case (boundary)
			case (.true.) ! closed
				call secant(p,bp)
			case (.false.) ! open
				call secant(p,bp)
		end select
  else
    call update_w(p,bp)
  end if
  steptime = steptime + 1 
  time = dt

!   !------------------ FOLLOWING TEMPORAL STEP ---------------------
  DO WHILE (time.LE. tfin)

		select case (boundary)
			case (.true.) ! closed
				call close_boundary(p,bp)
				call quantitiesatgridnodes_CLOSE(p,bp,gr)
				call gradient_correction_CLOSE(p,bp)
			case (.false.) ! open
			 	call quantitiesatgridnodes_OPEN(p,bp,gr)
			 	call gradient_correction_OPEN(p,bp)
		end select

    call locpar(p,bp,gr)
    call derivatives(gr)
    call pressure(gr)
    call avepressure(gr)
    call interp(p,bp,gr)

    if (boundary.eqv..true.) then ! closed
  		call boundarycupd(p,bp)
 	end if

    call turbulent_viscosity(p,gr)
    if (cfl_fixed.eq.0) then
        call stability(p,gr)
    else
        call CFL_routine(p,gr)
    endif

    call acceleration(p,bp)
    call leapfrog(p,bp)

	select case (boundary)
	case (.true.) !closed
			call collision(p)
		case (.false.) !open
			call change_particles(p,bp,gr)
            call open_boundary(p,bp,gr)
	end select

    if (SLupd.eqv..true.) then
        select case (boundary)
            case (.true.) ! closed
                call close_boundary(p,bp)!boundarycupd(p,bp)
                call secant(p,bp)
            case (.false.) ! open
                call secant(p,bp)
            end select
    else
        call update_w(p,bp)
    end if
    time = time +dt
    steptime = steptime + 1
    if(mod(steptime,res_freq).eq.0) then
      call output(gr,p,bp)
      print*,"tempo", time
    end if

  END DO

  CALL cpu_time(finish)                                                   ! final cpu time
  PRINT*,"CPU time:",finish-start
  OPEN(unit=100,file="test.dat",status = "old", position="append")
  write(100,*) "cpu time", finish-start
  close(100)
end program solver
