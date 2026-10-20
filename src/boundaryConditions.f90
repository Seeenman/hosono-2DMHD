module boundaryConditions


    use definitions, only: ndim, xdir, ydir, nfaces
    use simulation, only: sim_BC
    use gridBlock, only: gridBlock_t

    implicit none

    private

    public :: boundaryConditions_apply, boundaryConditions_applyBlk, boundaryConditions_applyMOOD

contains

    subroutine boundaryConditions_applyMOOD(blk)
        ! purpose:      Apply boundary conditions to the MOOD state of a block
        !               (blk%scheme and blk%recompute) so that guard cells
        !               match the cells they are images of. The detection step
        !               only lowers the order of interior cells, but the fluxes
        !               on the faces at the domain boundary, and the constrained
        !               transport corrections next to it, are built from
        !               reconstructions in the guard cells. If a guard cell is
        !               not lowered along with its image, the same physical face
        !               gets different fluxes on the two sides of a periodic
        !               boundary, which breaks conservation and div(B)=0 there.
        !
        ! Inputs:       - blk (gridBlock_t) the block whose MOOD state to fill
        !
        ! Outputs:      - blk%scheme and blk%recompute with guard cells and
        !                 their interior images made consistent
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk

        if (sim_BC == 'outflow') then
            ! guard cells are copies of the boundary cells rather than images
            ! of other interior cells, so there is nothing to keep consistent
            return
        else if (sim_BC == 'periodic') then
            call periodicMOOD(blk)
        else
            write(*,*) "=========================================================================="
            write(*,*) "Unrecognized choice of boundary condition: ", trim(sim_BC)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "sim_BC is set to a valid string, e.g., 'outflow' or 'periodic'"
            write(*,*) "=========================================================================="
            stop
        end if

    end subroutine boundaryConditions_applyMOOD

    subroutine boundaryConditions_applyBlk(blk, UorV)
        ! purpose:      Apply boundary conditions to the conservative or
        !               primitive variables living on a block. Convenience
        !               wrapper around boundaryConditions_apply that pulls
        !               the index bookkeeping out of the block and the
        !               number of variables out of the array itself.
        ! 
        ! Inputs:       - blk (gridBlock_t) the block the variables live on
        !               - UorV (real array) all conservative or primitive
        !                 variables at every cell of blk, e.g. grid_state%U,
        !                 grid_state%V, or the U of an RK substage
        !               
        ! Outputs:      - UorV (real array) with its guard cells filled
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in) :: blk
        ! UorV is only handed on to boundaryConditions_apply (which declares
        ! its bounds explicitly), so its bounds are not needed here.
        ! contiguous guarantees no temporary copy is made on that call.
        real, contiguous, intent(in out) :: UorV(:,:,:)

        call boundaryConditions_apply(UorV, SIZE(UorV, 1), blk%minIdx, blk%maxIdx, &
                                      blk%strtIdx, blk%stopIdx, blk%NGC)

    end subroutine boundaryConditions_applyBlk

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

    pure subroutine outflow(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)

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

    pure subroutine periodic(UorV, nVars, minIdx, maxIdx, strtIdx, stopIdx, NGC)
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

    pure subroutine periodicMOOD(blk)
        ! Give every guard cell and the interior cell it is a periodic image
        ! of the lower of their two orders, and flag both for recompute if
        ! either one is flagged. Flags have to go both ways: lowering an
        ! interior cell next to the boundary flags its neighbor across the
        ! boundary, which is a guard cell, so the interior image of that
        ! guard cell has to be recomputed as well.
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        ! local variables
        integer :: i, j, ip, jp
        logical :: eitherRecomp
        integer :: lowerScheme(nfaces)

        ! x direction. Every row is swept, guard rows included, so that the
        ! corner guard cells are reached when the y direction is swept below.
        do j=blk%minIdx(ydir), blk%maxIdx(ydir)
            do i=blk%minIdx(xdir), blk%maxIdx(xdir)
                if (i < blk%strtIdx(xdir)) then
                    ip = i + blk%N(xdir)
                else if (i > blk%stopIdx(xdir)) then
                    ip = i - blk%N(xdir)
                else
                    CYCLE ! interior cell, handled through its guard cell image
                end if
                ! propagate recomputeFLICD flag
                eitherRecomp = blk%recomputeFLICD(i,j) .or. blk%recomputeFLICD(ip,j)
                blk%recomputeFLICD(i, j) = eitherRecomp
                blk%recomputeFLICD(ip,j) = eitherRecomp
                ! propagate recompute flag
                eitherRecomp = blk%recompute(i,j) .or. blk%recompute(ip,j)
                blk%recompute(i, j) = eitherRecomp
                blk%recompute(ip,j) = eitherRecomp
                ! propagate scheme
                lowerScheme = MIN(blk%scheme(:,i,j), blk%scheme(:,ip,j))
                blk%scheme(:,i, j) = lowerScheme
                blk%scheme(:,ip,j) = lowerScheme
            end do
        end do

        ! y direction
        do j=blk%minIdx(ydir), blk%maxIdx(ydir)
            if (j < blk%strtIdx(ydir)) then
                jp = j + blk%N(ydir)
            else if (j > blk%stopIdx(ydir)) then
                jp = j - blk%N(ydir)
            else
                CYCLE ! interior row, handled through its guard row image
            end if
            do i=blk%minIdx(xdir), blk%maxIdx(xdir)
                ! propagate recomputeFLICD flag
                eitherRecomp = blk%recomputeFLICD(i,j) .or. blk%recomputeFLICD(i,jp)
                blk%recomputeFLICD(i,j ) = eitherRecomp
                blk%recomputeFLICD(i,jp) = eitherRecomp
                ! propagate recompute flag
                eitherRecomp = blk%recompute(i,j) .or. blk%recompute(i,jp)
                blk%recompute(i,j ) = eitherRecomp
                blk%recompute(i,jp) = eitherRecomp
                ! propagate scheme 
                lowerScheme = MIN(blk%scheme(:,i,j), blk%scheme(:,i,jp))
                blk%scheme(:,i,j ) = lowerScheme
                blk%scheme(:,i,jp) = lowerScheme
            end do
        end do

    end subroutine periodicMOOD

end module boundaryConditions
