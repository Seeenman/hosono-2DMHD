! Module with subroutines for setting the initial conditions
! for different simulations

module initialCondition

    use definitions
    use readParamFile, only: readParamFile_real, readParamFile_char
    use eos, only: eos_eintIdealGas
    use simulation, only: sim_gamma
    use boundaryConditions, only: boundaryConditions_applyBlk
    use convert, only: convert_prim2cons
    use gridBlock, only: gridBlock_t
    use gridState, only: gridState_t

    implicit none

    private

    public :: initialCondition_set

contains

    subroutine initialCondition_set(paramfile, blk, state)
        ! purpose:      Set the initial conditions
        !               
        ! 
        ! Inputs:       - paramfile (character) name of the parameter file
        !               - blk (gridBlock_t) the block to set the initial
        !                 condition on, carrying its own indices and geometry
        !               - state (gridState_t) the fluid state on blk
        !               
        ! Outputs:      - state%V (real array) primitive variables at every cell
        !               - state%U (real array) conservative variables at every cell
        ! ------------------------------------------------------------
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in) :: blk
        type(gridState_t), intent(in out) :: state

        character(len=max_string_length) :: IC_type
        integer i,j

        write(*,*) "=============================================================="
        write(*,*) "Setting initial conditions"
        write(*,*) "--------------------------------------------------------------"

        ! read the initial condition to be used from the parameter file
        IC_type = readParamFile_char(paramfile, "IC_type")
        ! set conservative variables
        if (trim(IC_type)=="explosion2d") then
            call explosion2d(paramfile, blk, state)
        else if (trim(IC_type)=="sedov2d") then
            call sedov2d(paramfile, blk, state)
        else if (trim(IC_type)=="OrszagTang2d") then
            call OrszagTang2D(paramfile, blk, state)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of initial condition: ", trim(IC_type)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "IC_type is set to a valid string, e.g., 'OrszagTang2D'"
            write(*,*) "=========================================================================="
            stop
        end if

        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                ! set internal energy and adiabatic index
                state%V(eint_var, i,j) = eos_eintIdealGas(state%V(pres_var,i,j), state%V(dens_var,i,j), sim_gamma)
                state%V(gamm_var, i,j) = sim_gamma

                ! also initialize conservative variables
                state%U(:,i,j) = convert_prim2cons(state%V(:,i,j))
            end do
        end do

        ! update boundary conditions because we have only set interior cells so far
        call boundaryConditions_applyBlk(blk, state%U)
        call boundaryConditions_applyBlk(blk, state%V)

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "Initial conditions set"
        write(*,*) "=============================================================="

    end subroutine initialCondition_set

    subroutine explosion2d(paramfile, blk, state)
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in) :: blk
        type(gridState_t), intent(in out) :: state

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

        state%V = 0.0
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                if (sqrt((blk%x(i)-shockCenter_x)**2 + (blk%y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    state%V(dens_var,i,j) = densIns
                    state%V(velx_var,i,j) = velxIns
                    state%V(vely_var,i,j) = velyIns
                    state%V(velz_var,i,j) = velzIns
                    state%V(pres_var,i,j) = presIns
                else
                    ! outside the circle
                    state%V(dens_var,i,j) = densOut
                    state%V(velx_var,i,j) = velxOut
                    state%V(vely_var,i,j) = velyOut
                    state%V(velz_var,i,j) = velzOut
                    state%V(pres_var,i,j) = presOut
                end if
            end do
        end do

    end subroutine explosion2d

    subroutine sedov2d(paramfile, blk, state)
        ! Initial conditions for 2D sedov test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in) :: blk
        type(gridState_t), intent(in out) :: state

        real :: shockCenter_x, shockCenter_y, shockRad
        real, dimension(nPrimVars) :: V_ins, V_out
        integer :: i,j, nn

        ! read required values from parameter file
        shockCenter_x = readParamFile_real(paramfile, "IC_shockCenter_x")
        shockCenter_y = readParamFile_real(paramfile, "IC_shockCenter_y")

        shockRad   = 3.5*MIN(blk%dl(xdir), blk%dl(ydir))

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

        state%V = 0.0
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                if (sqrt((blk%x(i)-shockCenter_x)**2 + (blk%y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    state%V(:, i,j) = V_ins
                else
                    ! outside the circle
                    state%V(:, i,j) = V_out
                end if
            end do
        end do

    end subroutine sedov2d

    subroutine OrszagTang2D(paramfile, blk, state)
        ! Initial conditions for 2D Orszag Tang mhd test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in) :: blk
        type(gridState_t), intent(in out) :: state

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
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                xx = blk%x(i)
                yy = blk%y(j)
                state%V(dens_var,i,j) =  dens
                state%V(velx_var,i,j) = -U0*SIN(pi*yy*2.0)
                state%V(vely_var,i,j) =  U0*SIN(pi*xx*2.0)
                state%V(velz_var,i,j) =  0.0
                state%V(magx_var,i,j) = -B0*SIN(pi*yy*2.0)
                state%V(magy_var,i,j) =  B0*SIN(pi*xx*4.0)
                state%V(magz_var,i,j) =  0.0
                state%V(pres_var,i,j) =  pres
            end do
        end do


    end subroutine OrszagTang2D

end module initialCondition
