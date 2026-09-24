module gridBlock

    ! Defines gridBlock_t: a collection of everything that describes one
    ! rectangular block of the grid, plus the fluid state living on it.
    !
    ! Why this type exists:
    !
    !   Routines that operate on a block (output, initial conditions,
    !   boundary conditions, the solver) otherwise need many separate 
    !   grid related arguments passed to them. Condensing them into gridBlock_t
    !   keeps subroutnine argument lists much shorter while still passing the grid data in
    !   explicitly.
    !
    ! What does NOT belong here:
    !
    !   Anything that is a property of the numerical method rather than of
    !   this piece of the domain - the GP stencil radius, the quadrature
    !   rule, the adiabatic index. Those are identical on every rank and
    !   live in the grid and simulation modules.

    use definitions, only: ndim, nConsVars, nPrimVars, xdir, ydir, nfaces

    implicit none

    private

    public :: gridBlock_t, gridBlock_alloc, gridBlock_dealloc

    type :: gridBlock_t

        ! ---- index bookkeeping ----
        integer :: NGC = 0                                   ! guard cells padding each side
        integer, dimension(ndim) :: N = 0                    ! number of interior cells
        integer, dimension(ndim) :: minIdx = 0, maxIdx = 0   ! index of first/last guard cell
        integer, dimension(ndim) :: strtIdx = 0, stopIdx = 0 ! index of first/last interior cell

        ! ---- geometry ----
        real, dimension(ndim) :: domainBeg = 0.0, domainEnd = 0.0 ! bounds of this block
        real, dimension(ndim) :: dl = 0.0                         ! holds dx and dy
        real, allocatable :: x(:), y(:)                           ! cell-center coordinates

        ! ---- fluid state ----
        ! dimensions correspond to (variable, xcoordinate, ycoordinate)
        real, allocatable :: U(:,:,:) ! conservative variables
        real, allocatable :: V(:,:,:) ! primitive variables

        ! ---- Riemann problem arrays ----
        ! Pointwise conservative variables reconstructed at face quadrature
        ! points, and the fluxes from the local Riemann problems there.
        ! (variable, direction, quadrature_point, xcoordinate, ycoordinate)
        real, allocatable :: faceVals(:,:,:,:,:)
        real, allocatable :: flux(:,:,:,:,:)

        ! ---- MOOD scheme cascade ----
        ! which reconstruction scheme to use at each face of each cell
        ! 3->7th order GP-R3
        ! 2->5th order GP-R2
        ! 1->3rd order GP-R1
        ! 0->1st order Godunov (FOG)
        ! (faceIdx, i, j)
        integer, allocatable :: scheme(:,:,:)

    end type gridBlock_t

contains

    pure subroutine gridBlock_alloc(blk, nquad)
        ! purpose:      Allocate and zero the arrays belonging to a block.
        !               The index bookkeeping and geometry components
        !               (N, minIdx, maxIdx, dl, ...) must already be set.
        !
        ! Inputs:       - blk (gridBlock_t) the block to allocate
        !               - nquad (integer) number of face quadrature points,
        !                 a property of the numerical method rather than of
        !                 the block, so it is passed in
        !
        ! Outputs:      - blk (gridBlock_t) with its arrays allocated and zeroed
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        integer, intent(in) :: nquad

        ! cell-center coordinates
        allocate(blk%x(blk%minIdx(xdir):blk%maxIdx(xdir)))
        allocate(blk%y(blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! conservative and primitive variables
        allocate(blk%U(nConsVars, &
                       blk%minIdx(xdir):blk%maxIdx(xdir), &
                       blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(blk%V(nPrimVars, &
                       blk%minIdx(xdir):blk%maxIdx(xdir), &
                       blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! face reconstructions and fluxes
        allocate(blk%faceVals(nConsVars, nquad, nfaces, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(blk%flux(nConsVars, ndim, nquad, &
                          blk%minIdx(xdir):blk%maxIdx(xdir), &
                          blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(blk%scheme(nfaces, &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))

        !!! zero all of these out !!!
        blk%x = 0.0
        blk%y = 0.0
        blk%U = 0.0
        blk%V = 0.0
        blk%faceVals = 0.0
        blk%flux = 0.0
        blk%scheme = 0

    end subroutine gridBlock_alloc

    pure subroutine gridBlock_dealloc(blk)
        ! purpose:      Deallocate the arrays belonging to a block
        !
        ! Inputs:       - blk (gridBlock_t) the block to deallocate
        !
        ! Outputs:      - none
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk

        deallocate(blk%x)
        deallocate(blk%y)
        deallocate(blk%U)
        deallocate(blk%V)
        deallocate(blk%faceVals)
        deallocate(blk%flux)
        deallocate(blk%scheme)

    end subroutine gridBlock_dealloc

end module gridBlock
