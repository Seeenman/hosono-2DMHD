module definitions

    implicit none

    ! grid values
    integer, parameter :: ndim = 2
    integer, parameter :: xdir = 1
    integer, parameter :: ydir = 2

    ! conservative variable indices
    integer, parameter :: nConsVars = 8 ! number of conservative variables
    integer, parameter :: dens_var = 1 ! density
    integer, parameter :: momx_var = 2 ! x-momentum
    integer, parameter :: momy_var = 3 ! y-momentum
    integer, parameter :: momz_var = 4 ! z-momentum
    integer, parameter :: ener_var = 5 ! total energy
    integer, parameter :: magx_var = 6 ! x-magnetic field
    integer, parameter :: magy_var = 7 ! y-magnetic field
    integer, parameter :: magz_var = 8 ! z-magnetic field

    ! primitive variable indices
    integer, parameter :: nPrimVars = 10 ! number of conservative variables
    integer, parameter :: velx_var = 2 ! x-velocity
    integer, parameter :: vely_var = 3 ! y-velocity
    integer, parameter :: velz_var = 4 ! z-velocity
    integer, parameter :: pres_var = 5 ! pressure
    integer, parameter :: eint_var = 9 ! internal energy
    integer, parameter :: gamm_var = 10 ! adiabatic index

    ! other
    integer, parameter :: max_string_length = 800
    real, parameter :: pi = ACOS(-1.0)

contains

end module definitions
