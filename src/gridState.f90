module gridState

    ! Defines gridState_t: the fluid state living on one rectangular block
    ! of the grid (see gridBlock.f90 for the block itself).
    !
    ! Why this is separate from gridBlock_t:
    !
    !   A multi-stage time integrator (e.g. SSP-RK3) needs several copies
    !   of the conservative variables at once (U^n, U^(1), U^(2), ...), but
    !   only one copy of the block's geometry, index bookkeeping, and
    !   scratch arrays (faceVals, flux, scheme). Keeping the solution out of
    !   the block lets the routines used inside an RK stage (reconstruction,
    !   forward Euler, boundary conditions) take a plain U array, so they
    !   can be handed either grid_state%U or the U of a substage.
    !
    !   The substages themselves only need U, so they are plain allocatable
    !   arrays with the same bounds as gridState_t%U rather than full
    !   gridState_t copies. gridState_t is used for the solution of the run,
    !   where U and V travel together (output, cfl, initial conditions).
    !
    ! A gridState_t carries no index information of its own. Its arrays
    ! are dimensioned from the gridBlock_t it was allocated against, and
    ! it should always be passed around together with that block.

    use definitions, only: nConsVars, nPrimVars, xdir, ydir
    use gridBlock, only: gridBlock_t

    implicit none

    private

    public :: gridState_t, gridState_alloc, gridState_dealloc

    type :: gridState_t

        ! dimensions correspond to (variable, xcoordinate, ycoordinate)
        ! and include guard cells
        real, allocatable :: U(:,:,:) ! conservative variables
        real, allocatable :: V(:,:,:) ! primitive variables

    end type gridState_t

contains

    pure subroutine gridState_alloc(state, blk)
        ! purpose:      Allocate and zero the arrays belonging to a state.
        !
        ! Inputs:       - state (gridState_t) the state to allocate
        !               - blk (gridBlock_t) the block the state lives on.
        !                 Its index bookkeeping (minIdx, maxIdx) must
        !                 already be set.
        !
        ! Outputs:      - state (gridState_t) with its arrays allocated and zeroed
        ! ------------------------------------------------------------
        implicit none
        type(gridState_t), intent(in out) :: state
        type(gridBlock_t), intent(in) :: blk

        ! conservative and primitive variables
        allocate(state%U(nConsVars, &
                         blk%minIdx(xdir):blk%maxIdx(xdir), &
                         blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(state%V(nPrimVars, &
                         blk%minIdx(xdir):blk%maxIdx(xdir), &
                         blk%minIdx(ydir):blk%maxIdx(ydir)))

        !!! zero all of these out !!!
        state%U = 0.0
        state%V = 0.0

    end subroutine gridState_alloc

    pure subroutine gridState_dealloc(state)
        ! purpose:      Deallocate the arrays belonging to a state
        !
        ! Inputs:       - state (gridState_t) the state to deallocate
        !
        ! Outputs:      - none
        ! ------------------------------------------------------------
        implicit none
        type(gridState_t), intent(in out) :: state

        deallocate(state%U)
        deallocate(state%V)

    end subroutine gridState_dealloc

end module gridState
