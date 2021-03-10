module global_var

! -------------------------------
! -------------------------------
! Multilayer shallow water solver
! by Luca Macavero
! -------------------------------
! -------------------------------

  ! constants & parameters
  ! -pi: pi greco
  ! -g: gravitational field [N/kg]
  ! -dt: time step [s]
  ! -tfin: final simulation time [s]
  ! -htot: initial height of the film [m]
  ! -patm: atmospheric pressure [Pa]
  ! -nu,sigma,rho: properties of the liquid [m^2/s,N/m,kg/m^3]
  ! -beta: angle of inclination of the wall [rad]
  ! -domain: array of the extremes of the domain [m]
  ! -nlayer: number of layers
  ! -numpar: initial number of particles
  ! -steptime: integer timestepimplicit none
  ! -ninterp: number of interpolation nodes
  ! -hsx: initial height at the left boundary
  ! -hdx: initial height at the right boundary
  real*8,parameter :: pi = 4.*atan(1.)
  real*8,parameter :: g = 9.80665
  real*8 :: dt
  real*8 :: tfin
  real*8 :: htot
  real*8 :: patm
  real*8 :: nu,sigma,rho
  real*8 :: beta
  real*8 :: hsx,hdx,inVsx,inVdx
  real*8 :: domain(2)
  real*8,allocatable :: inusx(:),inudx(:)
  integer :: idx_end_dom,idx_init_dom

  integer :: numpar
  integer :: nlayer
  integer :: steptime
  integer :: ninterp

  integer :: nbp
  integer :: testcase
  logical :: boundary, SLupd, turbulence
  integer :: res_freq
  real*8 :: art_visc_coeff  
  real*8 :: cfl_fixed 
  
  integer :: np ! serve solo nel caso di close bc

  ! contains particles information
  type particles
    real*8,allocatable :: x(:),u(:),a(:), nu_t(:)
    real*8,allocatable :: V(:),w(:),have(:),l(:),B(:),inhave(:),inl(:)
    real*8,allocatable :: xold(:),uold(:),aold(:)

    real*8,allocatable :: pave(:),pplus(:),pmin(:),dudzplus(:),dudzmin(:)
    real*8,allocatable :: dhdxplus(:),dhdxmin(:)

    integer,allocatable :: loc(:)
    integer :: np
  end type particles

  type mesh
    real*8,allocatable :: x(:),u(:),have(:),dudz(:),dhdx(:),l(:),pint(:)
    real*8,allocatable :: pave(:),zcoord(:)
  end type mesh

contains

  subroutine read_input()
    implicit none
    real*8 :: length


    namelist /probl_domain/testcase,beta,length,boundary,nlayer,tfin
    namelist /SPH/numpar,ninterp,dt,cfl_fixed,art_visc_coeff,SLupd
    namelist /physics/patm,rho,nu, turbulence
    namelist /results/res_freq
    open(0, file="file_input.txt")
    read(0, nml= probl_domain)
    read(0, nml = SPH)
    read(0, nml = physics)
    read(0, nml = results)
    close(0)
    domain(1) = 0
    domain(2) = length
    beta = beta*pi/180.
  end subroutine read_input
end module global_var
