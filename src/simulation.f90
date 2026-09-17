module simulation

    use readParamFile, only: readParamFile_int, readParamFile_real, readParamFile_char, readParamFile_logical
    use definitions, only: max_string_length

    implicit none

    ! eos
    real :: sim_gamma

    ! floor values
    real :: sim_smallB, sim_smallDensity, sim_smallEnergy, sim_smallPressure

    ! solver options
    character(len=max_string_length) :: sim_riemannSolver, sim_BC, &
                                        sim_dataFileBaseName

    ! output type
    logical :: sim_outputHdf5, sim_outputAscii, sim_forceHydro

    ! maximum number of time steps, output frequency in time step
    integer :: sim_nstepmax, sim_outputFreqStep

    ! courant number, tmax, output frequency in time
    real :: sim_cfl, sim_tmax, &
                      sim_outputFreqTime

contains

    subroutine simulation_init(paramfile)
        implicit none
        character(len=max_string_length), intent(in) :: paramfile

        write(*,*) "=============================================================="
        write(*,*) "Initializing simulation"
        write(*,*) "--------------------------------------------------------------"

        sim_smallB = readParamFile_real(paramfile, 'sim_smallB')
        sim_smallPressure = readParamFile_real(paramfile, 'sim_smallPressure')
        sim_smallDensity = readParamFile_real(paramfile, 'sim_smallDensity')
        sim_smallEnergy = readParamFile_real(paramfile, 'sim_smallEnergy')
        sim_gamma = readParamFile_real(paramfile, 'sim_gamma')
        sim_riemannSolver = readParamFile_char(paramfile, 'sim_riemannSolver')
        sim_cfl = readParamFile_real(paramfile, 'sim_cfl')
        sim_tmax = readParamFile_real(paramfile, 'sim_tmax')
        sim_nstepmax = readParamFile_int(paramfile, 'sim_nstepmax')
        sim_BC = readParamFile_char(paramfile, 'sim_BC')
        sim_dataFileBaseName = readParamFile_char(paramfile, 'sim_dataFileBaseName')
        sim_outputFreqTime = readParamFile_real(paramfile, 'sim_outputFreqTime')
        sim_outputFreqStep = readParamFile_int(paramfile, 'sim_outputFreqStep')
        sim_outputHdf5 = readParamFile_logical(paramfile, 'sim_outputHdf5')
        sim_outputAscii = readParamFile_logical(paramfile, 'sim_outputAscii')
        sim_forceHydro = readParamFile_logical(paramfile, 'sim_forceHydro')

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "Simulation initialized"
        write(*,*) "=============================================================="

    end subroutine simulation_init

    subroutine simulation_finalize()
        implicit none
        write(*,*) "=============================================================="
        write(*,*) "Simulation variables deallocated (nothing to deallocate)."
        write(*,*) "=============================================================="
    end subroutine simulation_finalize


end module simulation
