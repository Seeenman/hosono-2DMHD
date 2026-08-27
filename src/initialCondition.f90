! Module with subroutines for setting the initial conditions
! for different simulations

module initialCondition

    use definitions
    use readParamFile, only: readParamFile_real, readParamFile_char
    use eos, only: eos_eintIdealGas
    use simulation, only: sim_gamma
    use boundaryConditions, only: boundaryConditions_apply
    use convert, only: convert_prim2cons

    implicit none

    private

    public :: initialCondition_set

contains

    subroutine initialCondition_set(paramfile, U, V, N, minIdx, maxIdx,&
                                    strtIdx, stopIdx, NGC, dl, x, y)
        ! purpose:      Set the initial conditions
        !               
        ! 
        ! Inputs:       - U (real array) all conservative variables at every cell
        !               - V (real array) all primitive variables at every cell
        !               - N (integer array) number of cells in each direction
        !               - minIdx/maxIdx (integer arrays) index of first/last
        !                 guard cells in each direction
        !               - strtIdx/stopIdx (integer arrays) index of first/last
        !                 interior cell in each direction
        !               - NGC (integer) number of guard cells in each direction
        !               - dl (real array) holds dx and dy
        !               - x/y (real arrays) coordinates of cell centers
        !               
        ! Outputs:      - U (real array) all conservative variables at every cell
        ! ------------------------------------------------------------
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        integer, dimension(ndim), intent(in) :: minIdx, maxIdx, strtIdx, stopIdx
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC), dl(ndim)
        real, intent(in out) :: U(nConsVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))
        real, intent(in out) :: V(nPrimVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        character(len=max_string_length) :: IC_type
        integer i,j

        write(*,*) "=============================================================="
        write(*,*) "Setting initical conditions"
        write(*,*) "--------------------------------------------------------------"

        ! read the initial condition to be used from the parameter file
        IC_type = readParamFile_char(paramfile, "IC_type")
        ! set conservative variables
        if (trim(IC_type)=="explosion2d") then
            call explosion2d(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        else if (trim(IC_type)=="sedov2d") then
            call sedov2d(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, dl, x, y)
        else if (trim(IC_type)=="OrszagTang2d") then
            call OrszagTang2D(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of initial condition: ", trim(IC_type)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "IC_type is set to a valid string, e.g., 'OrszagTang2D'"
            write(*,*) "=========================================================================="
            stop
        end if

        ! set internal energy
        do j=strtIdx(YDIR), stopIdx(YDIR)
            do i=strtIdx(XDIR), stopIdx(XDIR)
                V(eint_var, i,j) = eos_eintIdealGas(V(pres_var,i,j), V(dens_var,i,j), sim_gamma)
            end do
        end do

        ! update boundary conditions because we have only set interior cells so far
        call boundaryConditions_apply(V, nPrimVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)

        ! also initialize conservative variables.
        ! note loops go over all cells including guard cells
        do j=minIdx(ydir), maxIdx(ydir)
            do i=minIdx(xdir), maxIdx(xdir)
                U(:,i,j) = convert_prim2cons(V(:,i,j))
            end do
        end do

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "Initial conditions set"
        write(*,*) "=============================================================="

    end subroutine initialCondition_set

    subroutine explosion2d(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC)
        real, intent(in out) :: V(nPrimVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        real :: shockCenter_x, shockCenter_y, shockRad, &
                    densIns, densOut, &
                    velxIns, velxOut, &
                    velyIns, velyOut, &
                    velzIns, velzOut, &
                    presIns, presOut
        integer :: i,j

        ! read required values from parameter file
        shockCenter_x = readParamFile_real(paramfile, "IC_shockCenter_x")
        shockCenter_y = readParamFile_real(paramfile, "IC_shockCenter_y")
        shockRad = readParamFile_real(paramfile, "IC_shockRad")
        densIns = readParamFile_real(paramfile, "IC_densIns")
        densOut = readParamFile_real(paramfile, "IC_densOut")
        velxIns = readParamFile_real(paramfile, "IC_velxIns")
        velyIns = readParamFile_real(paramfile, "IC_velyIns")
        velzIns = readParamFile_real(paramfile, "IC_velzIns")
        velxOut = readParamFile_real(paramfile, "IC_velxOut")
        velyOut = readParamFile_real(paramfile, "IC_velyOut")
        velzOut = readParamFile_real(paramfile, "IC_velzOut")
        presIns = readParamFile_real(paramfile, "IC_presIns")
        presOut = readParamFile_real(paramfile, "IC_presOut")

        V = 0.0
        do j=strtIdx(ydir), stopIdx(ydir)
            do i=strtIdx(xdir), stopIdx(xdir)
                if (sqrt((x(i)-shockCenter_x)**2 + (y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    V(dens_var,i,j) = densIns
                    V(velx_var,i,j) = velxIns
                    V(vely_var,i,j) = velyIns
                    V(velz_var,i,j) = velzIns
                    V(ener_var,i,j) = presIns
                else
                    ! outside the circle
                    V(dens_var,i,j) = densOut
                    V(momx_var,i,j) = velxOut*densOut
                    V(momy_var,i,j) = velyOut*densOut
                    V(momz_var,i,j) = velzOut*densOut
                    V(ener_var,i,j) = presOut
                end if
            end do
        end do

    end subroutine explosion2d

    subroutine sedov2d(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, dl, x, y)
        ! Initial conditions for 2D sedov test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC), dl(ndim)
        real, intent(in out) :: V(nPrimVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        real :: shockCenter_x, shockCenter_y, shockRad
        real, dimension(nPrimVars) :: V_ins, V_out
        integer :: i,j, nn

        ! read required values from parameter file
        shockCenter_x = readParamFile_real(paramfile, "IC_shockCenter_x")
        shockCenter_y = readParamFile_real(paramfile, "IC_shockCenter_y")

        shockRad   = 3.5*MIN(dl(xdir), dl(ydir))

        ! cylindricaly geometry
        nn = 2

        ! primitive variables inside of circle
        V_ins(dens_var) = 1.0
        V_ins(velx_var) = 0.0
        V_ins(vely_var) = 0.0
        V_ins(velz_var) = 0.0
        V_ins(pres_var) = 3*(sim_gamma-1)*1.0/((nn+1)*pi*shockRad**nn)
        ! primitive variables outside of circle
        V_out(dens_var) = 1.0
        V_out(velx_var) = 0.0
        V_out(vely_var) = 0.0
        V_out(velz_var) = 0.0
        V_out(pres_var) = 1e-5

        V = 0.0
        do j=strtIdx(ydir), stopIdx(ydir)
            do i=strtIdx(xdir), stopIdx(xdir)
                if (sqrt((x(i)-shockCenter_x)**2 + (y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    V(:, i,j) = V_ins
                else
                    ! outside the circle
                    V(:, i,j) = V_out
                end if
            end do
        end do

    end subroutine sedov2d

    subroutine OrszagTang2D(paramfile, V, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        ! Initial conditions for 2D Orszag Tang mhd test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC)
        real, intent(in out) :: V(nPrimVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        real :: U0, dens
        real :: B0, pres
        integer :: i,j
        real :: xx, yy

        ! read required values from parameter file
        u0   = readParamFile_real(paramfile, "IC_U0")
        dens = readParamFile_real(paramfile, "IC_dens")

        ! set other values accordingly
        B0   = 1/sim_gamma
        pres = 1/sim_gamma

        ! Set primitive variables
        do j=strtIdx(ydir), stopIdx(ydir)
            do i=strtIdx(xdir), stopIdx(xdir)
                xx = x(i)
                yy = y(j)
                V(dens_var,i,j) =  dens
                V(velx_var,i,j) = -U0*SIN(pi*yy*2.0)
                V(vely_var,i,j) =  U0*SIN(pi*xx*2.0)
                V(velz_var,i,j) =  0.0
                V(magx_var,i,j) = -B0*SIN(pi*yy*2.0)
                V(magy_var,i,j) =  B0*SIN(pi*xx*4.0)
                V(magz_var,i,j) =  0.0
                V(pres_var,i,j) =  pres
            end do
        end do


    end subroutine OrszagTang2D

end module initialCondition
