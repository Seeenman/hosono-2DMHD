! Module with subroutines for setting the initial conditions
! for different simulations

module initialCondition

    use definitions
    use readParamFile, only: readParamFile_real, readParamFile_char
    use eos, only: eos_eintIdealGas
    use simulation, only: sim_gamma
    use boundaryConditions, only: boundaryConditions_apply
    use convert, only: convert_prim2cons
    use gridBlock, only: gridBlock_t

    implicit none

    private

    public :: initialCondition_set

contains

    subroutine initialCondition_set(paramfile, blk)
        ! purpose:      Set the initial conditions
        !               
        ! 
        ! Inputs:       - paramfile (character) name of the parameter file
        !               - blk (gridBlock_t) the block to set the initial
        !                 condition on, carrying its own indices, geometry,
        !                 and fluid state
        !               
        ! Outputs:      - blk%V (real array) primitive variables at every cell
        !               - blk%U (real array) conservative variables at every cell
        ! ------------------------------------------------------------
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in out) :: blk

        character(len=max_string_length) :: IC_type
        integer i,j

        write(*,*) "=============================================================="
        write(*,*) "Setting initial conditions"
        write(*,*) "--------------------------------------------------------------"

        ! read the initial condition to be used from the parameter file
        IC_type = readParamFile_char(paramfile, "IC_type")
        ! set conservative variables
        if (trim(IC_type)=="explosion2d") then
            call explosion2d(paramfile, blk)
        else if (trim(IC_type)=="sedov2d") then
            call sedov2d(paramfile, blk)
        else if (trim(IC_type)=="OrszagTang2d") then
            call OrszagTang2D(paramfile, blk)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of initial condition: ", trim(IC_type)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "IC_type is set to a valid string, e.g., 'OrszagTang2D'"
            write(*,*) "=========================================================================="
            stop
        end if

        ! set internal energy and adiabatic index
        do j=blk%strtIdx(YDIR), blk%stopIdx(YDIR)
            do i=blk%strtIdx(XDIR), blk%stopIdx(XDIR)
                blk%V(eint_var, i,j) = eos_eintIdealGas(blk%V(pres_var,i,j), blk%V(dens_var,i,j), sim_gamma)
                blk%V(gamm_var, i,j) = sim_gamma
            end do
        end do

        ! update boundary conditions because we have only set interior cells so far
        call boundaryConditions_apply(blk%V, nPrimVars, blk%minIdx, blk%maxIdx, &
                                      blk%strtIdx, blk%stopIdx, blk%NGC)

        ! also initialize conservative variables.
        ! note loops go over all cells including guard cells
        do j=blk%minIdx(ydir), blk%maxIdx(ydir)
            do i=blk%minIdx(xdir), blk%maxIdx(xdir)
                blk%U(:,i,j) = convert_prim2cons(blk%V(:,i,j))
            end do
        end do

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "Initial conditions set"
        write(*,*) "=============================================================="

    end subroutine initialCondition_set

    subroutine explosion2d(paramfile, blk)
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in out) :: blk

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

        blk%V = 0.0
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                if (sqrt((blk%x(i)-shockCenter_x)**2 + (blk%y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    blk%V(dens_var,i,j) = densIns
                    blk%V(velx_var,i,j) = velxIns
                    blk%V(vely_var,i,j) = velyIns
                    blk%V(velz_var,i,j) = velzIns
                    blk%V(pres_var,i,j) = presIns
                else
                    ! outside the circle
                    blk%V(dens_var,i,j) = densOut
                    blk%V(velx_var,i,j) = velxOut
                    blk%V(vely_var,i,j) = velyOut
                    blk%V(velz_var,i,j) = velzOut
                    blk%V(pres_var,i,j) = presOut
                end if
            end do
        end do

    end subroutine explosion2d

    subroutine sedov2d(paramfile, blk)
        ! Initial conditions for 2D sedov test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in out) :: blk

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

        blk%V = 0.0
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                if (sqrt((blk%x(i)-shockCenter_x)**2 + (blk%y(j)-shockCenter_y)**2) <= shockRad) then
                    ! inside the circle
                    blk%V(:, i,j) = V_ins
                else
                    ! outside the circle
                    blk%V(:, i,j) = V_out
                end if
            end do
        end do

    end subroutine sedov2d

    subroutine OrszagTang2D(paramfile, blk)
        ! Initial conditions for 2D Orszag Tang mhd test problem
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        type(gridBlock_t), intent(in out) :: blk

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
                blk%V(dens_var,i,j) =  dens
                blk%V(velx_var,i,j) = -U0*SIN(pi*yy*2.0)
                blk%V(vely_var,i,j) =  U0*SIN(pi*xx*2.0)
                blk%V(velz_var,i,j) =  0.0
                blk%V(magx_var,i,j) = -B0*SIN(pi*yy*2.0)
                blk%V(magy_var,i,j) =  B0*SIN(pi*xx*4.0)
                blk%V(magz_var,i,j) =  0.0
                blk%V(pres_var,i,j) =  pres
            end do
        end do


    end subroutine OrszagTang2D

end module initialCondition
