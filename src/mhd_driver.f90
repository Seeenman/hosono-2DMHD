program mhd_driver

    use, intrinsic :: ieee_arithmetic, only: IEEE_IS_NAN

    use definitions, only: max_string_length
    use simulation, only: simulation_init, simulation_finalize, sim_tmax, sim_nstepmax
    use GP, only: GP_init, GP_finalize, GP_nQuadrature, GP_radius
    use grid, only: grid_init, grid_finalize, grid_block
    use initialCondition, only: initialCondition_set
    use output, only: output_write
    use cfl, only: cfl_computedt

    implicit none

    character(len=max_string_length) :: paramfile

    ! time stepping / output bookkeeping
    integer :: nStep, lastOutputStep, outputCounter
    real :: t, dt, lastOutputTime

    ! get name of parameter file
    call GET_COMMAND_ARGUMENT(1,paramfile)
    if (len_trim(paramfile)==0) then
        write(*,*) "=============================================================="
        write(*,*) "No parameter file specified."
        write(*,*) "Please specify a parameter file, e.g."
        write(*,*) "./hosono-2DMHD.ex myfile.par"
        write(*,*) "=============================================================="
        stop
    end if
    write(*,*) "=============================================================="
    write(*,*) "Reading input parameters from ", trim(paramfile)
    write(*,*) "=============================================================="

    ! ----------------------------------------------------
    ! set the variables defined inside of the simulation
    ! module by reading them from the parameter file
    ! ----------------------------------------------------
    call simulation_init(paramfile)

    ! ----------------------------------------------------
    ! Set the variables defined inside of the GP
    ! module
    ! ----------------------------------------------------
    call GP_init(paramfile)

    ! ----------------------------------------------------
    ! Allocate the global grid variables and arrays
    ! in grid_block 
    ! ----------------------------------------------------
    call grid_init(paramfile, GP_radius, GP_nQuadrature)

    ! ----------------------------------
    ! set initial conditions
    ! ----------------------------------
    call initialCondition_set(paramfile, grid_block)

    ! ----------------------------------
    ! Advance the solution in time up
    ! to tmax or nstepmax
    ! ----------------------------------
    nStep = 0
    t = 0.0
    dt = cfl_computedt(grid_block)
    call validTimeStep(dt, t, sim_tmax)
    lastOutputStep = 0
    lastOutputTime = 0.0
    outputCounter = 0
    
    ! write initial conditions to disk
    call output_write(nStep, t, dt, lastOutputStep, lastOutputTime, outputCounter, .true., grid_block)

    ! the main loop
    ! do while ((t < sim_tmax) .and. (nStep < sim_nstepmax))
    !     ! update dt based on cfl
    !     dt = cfl_computedt(grid_block)
    !     call validTimeStep(dt, t, sim_tmax)

    ! end do

    ! ----------------------------------------------------
    ! finalize (deallocate data)
    ! ----------------------------------------------------
    call grid_finalize()
    call GP_finalize()
    call simulation_finalize()

    write(*,*) "=============================================================="
    write(*,*) "Simulation has ended."
    write(*,*) "=============================================================="

contains

    subroutine validTimeStep(delta_t, current_time, tmax)
        ! Checks
        ! 1. if t+dt is less than tmax. If not, set dt=tmax-t
        ! 2. if dt is positive or NAN. If so, quit simulation.
        implicit none
        real, intent(in out) :: delta_t
        real, intent(in) :: current_time, tmax

        if (t+delta_t > tmax) then
            delta_t = tmax - current_time
        else if (delta_t < 0.0) then
            write(*,*) "=========================================================================="
            write(*,*) "dt is negative. Aborting simulation."
            write(*,*) "=========================================================================="
            stop
        else if (IEEE_IS_NAN(delta_t)) then
            write(*,*) "=========================================================================="
            write(*,*) "dt is NAN. Aborting simulation."
            write(*,*) "=========================================================================="
            stop
        end if

    end subroutine validTimeStep
    
end program mhd_driver
