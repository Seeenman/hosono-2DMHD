module boundaryConditions


    use definitions, only: ndim, xdir, ydir
    use readParamFile, only: readParamFile_char
    use simulation, only: sim_BC

    implicit none

    private

    public :: boundaryConditions_apply

contains

    subroutine boundaryConditions_apply(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)
        ! purpose:      Apply boundary conditions to the 
        !               conservative or primitive variables
        ! 
        ! Inputs:       - UorV (real array) all conservative or primitive variables
        !                 at every cell
        !               - N (integer) nConsVars or nPrimVars
        !               - N (integer array) number of cells in each direction
        !               - minIdx/maxIdx (integer arrays) index of first/last
        !                 guard cells in each direction
        !               - strtIdx/stopIdx (integer arrays) index of first/last
        !                 interior cell in each direction
        !               - NGC (integer) number of guard cells in each direction
        !               
        ! Outputs:      - UorV (real array) all conservative or primnitive variables
        !                 at every cell
        ! ------------------------------------------------------------
        implicit none
        integer, dimension(ndim), intent(in) :: minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC, nVars
        real, intent(in out) :: UorV(nVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        if (sim_BC == 'outflow') then
            call outflow(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)
        else if (sim_BC == 'periodic') then
            call periodic(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of boundary condition: ", trim(sim_BC)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "sim_BC is set to a valid string, e.g., 'outflow' or 'periodic'"
            write(*,*) "=========================================================================="
            stop
        end if

    end subroutine boundaryConditions_apply

    subroutine outflow(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)

        implicit none

        integer, dimension(ndim), intent(in) :: minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC, nVars
        real, intent(in out) :: UorV(nVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))
        ! local variables
        integer :: i

        ! lower x boundary condition
        do i=minIdx(xdir),NGC
            UorV(:, i, :) = &
                UorV(:, strtIdx(xdir), :)
        end do
        ! upper x boundary condition
        do i=stopIdx(xdir)+1,maxIdx(xdir)
            UorV(:, i, :) = &
                UorV(:, stopIdx(xdir), :)
        end do

        ! lower y boundary condition
        do i=minIdx(ydir),NGC
            UorV(:, :, i) = &
                UorV(:, :, strtIdx(ydir))
        end do
        ! upper y boundary condition
        do i=stopIdx(ydir)+1,maxIdx(ydir)
            UorV(:, :, i) = &
                UorV(:, :, stopIdx(ydir))
        end do

    end subroutine outflow

    subroutine periodic(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)
        implicit none

        integer, dimension(ndim), intent(in) :: minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC, nVars
        real, intent(in out) :: UorV(nVars, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        integer :: i

        ! lower x boundary condition
        do i=1,NGC
            UorV(:, strtIdx(xdir)-i, :) = &
                UorV(:, stopIdx(xdir)-i+1, :)
        end do
        ! upper x boundary condition
        do i=1,NGC
            UorV(:, stopIdx(xdir)+i, :) = &
                UorV(:, strtIdx(xdir)+i-1, :)
        end do

        ! lower y boundary condition
        do i=1,NGC
            UorV(:, :, strtIdx(ydir)-i) = &
                UorV(:, :, stopIdx(ydir)-i+1)
        end do
        ! upper y boundary condition
        do i=1,NGC
            UorV(:, :, stopIdx(ydir)+i) = &
                UorV(:, :, strtIdx(ydir)+i-1)
        end do

    end subroutine periodic

end module boundaryConditions
