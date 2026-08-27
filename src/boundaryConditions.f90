module boundaryConditions


    use definitions, only: ndim, xdir, ydir
    use readParamFile, only: readParamFile_char
    use simulation, only: sim_BC

    implicit none

    private

    public :: boundaryConditions_apply

contains

    subroutine boundaryConditions_apply(UorV, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        ! purpose:      Apply boundary conditions to the 
        !               conservative or primitive variables
        ! 
        ! Inputs:       - UorV (real array) all conservative or primitive variables
        !                 at every cell
        !               - N (integer array) number of cells in each direction
        !               - minIdx/maxIdx (integer arrays) index of first/last
        !                 guard cells in each direction
        !               - strtIdx/stopIdx (integer arrays) index of first/last
        !                 interior cell in each direction
        !               - NGC (integer) number of guard cells in each direction
        !               - x/y (real arrays) coordinates of cell centers
        !               
        ! Outputs:      - UorV (real array) all conservative or primnitive variables
        !                 at every cell
        ! ------------------------------------------------------------
        implicit none
        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC)
        real, intent(in out) :: UorV(:, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        if (sim_BC == 'outflow') then
            call outflow(UorV, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        else if (sim_BC == 'periodic') then
            call periodic(UorV, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of boundary condition: ", trim(sim_BC)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "sim_BC is set to a valid string, e.g., 'outflow' or 'periodic'"
            write(*,*) "=========================================================================="
            stop
        end if

    end subroutine boundaryConditions_apply

    subroutine outflow(U, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)

        implicit none

        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC)
        real, intent(in out) :: U(:, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))
        ! local variables
        integer :: i

        ! lower x boundary condition
        do i=minIdx(xdir),NGC
            U(:, i, :) = &
                U(:, strtIdx(xdir), :)
        end do
        ! upper x boundary condition
        do i=stopIdx(xdir)+1,maxIdx(xdir)
            U(:, i, :) = &
                U(:, stopIdx(xdir), :)
        end do

        ! lower y boundary condition
        do i=minIdx(ydir),NGC
            U(:, :, i) = &
                U(:, :, strtIdx(ydir))
        end do
        ! upper y boundary condition
        do i=stopIdx(ydir)+1,maxIdx(ydir)
            U(:, :, i) = &
                U(:, :, stopIdx(ydir))
        end do

    end subroutine outflow

    subroutine periodic(UorV, N, minIdx, maxIdx, strtIdx, stopIdx, NGC, x, y)
        implicit none

        integer, dimension(ndim), intent(in) :: N, minIdx, maxIdx, strtIdx, stopIdx
        integer, intent(in) :: NGC
        real, intent(in) :: x(N(xdir)+2*NGC), y(N(ydir)+2*NGC)
        real, intent(in out) :: UorV(:, &
                                  minIdx(xdir):maxIdx(xdir), &
                                  minIdx(ydir):maxIdx(ydir))

        integer :: i

        ! lower x boundary condition
        do i=1,NGC
            U(:, strtIdx(xdir)-i, :) = &
                U(:, stopIdx(xdir)-i+1, :)
        end do
        ! upper x boundary condition
        do i=1,NGC
            U(:, stopIdx(xdir)+i, :) = &
                U(:, strtIdx(xdir)+i-1, :)
        end do

        ! lower y boundary condition
        do i=1,NGC
            U(:, :, strtIdx(ydir)-i) = &
                U(:, :, stopIdx(ydir)-i+1)
        end do
        ! upper y boundary condition
        do i=1,NGC
            U(:, :, stopIdx(ydir)+i) = &
                U(:, :, strtIdx(ydir)+i-1)
        end do

    end subroutine periodic

end module boundaryConditions
