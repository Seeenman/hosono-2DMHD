program mhd_driver

    use definitions, only: max_string_length
    use simulation, only: simulation_init
    use grid
    use initialCondition, only: initialCondition_set
    use output, only: output_write

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
    ! Set the variables defined inside of the grid
    ! module by reading them from the parameter file.
    ! Also allocate the global grid arrays.
    ! ----------------------------------------------------
    call grid_init(paramfile)

    ! ----------------------------------
    ! set initial conditions
    ! ----------------------------------
    call initialCondition_set(paramfile, grid_block)

    ! ----------------------------------
    ! write the initial condition to disk
    ! ----------------------------------
    nStep = 0
    t = 0.0
    dt = 0.0
    lastOutputStep = 0
    lastOutputTime = 0.0
    outputCounter = 0
    call output_write(nStep, t, dt, lastOutputStep, lastOutputTime, outputCounter, .true., grid_block)

    ! ----------------------------------------------------
    ! finalize (deallocate data)
    ! ----------------------------------------------------
    call grid_finalize()

    write(*,*) "=============================================================="
    write(*,*) "Simulation has ended."
    write(*,*) "=============================================================="
    
end program mhd_driver
