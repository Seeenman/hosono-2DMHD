module gridBlock

    ! Defines gridBlock_t: a collection of everything that describes one
    ! rectangular block of the grid, plus the scratch arrays used when
    ! advancing the fluid state living on it. The fluid state itself (U, V)
    ! is kept separately in gridState_t (see gridState.f90) so that a
    ! single block can be paired with several copies of U, e.g. the
    ! substages of a Runge-Kutta method.
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

    use definitions, only: ndim, nConsVars, xdir, ydir, nfaces

    implicit none

    private

    public :: gridBlock_t, gridBlock_alloc, gridBlock_dealloc, gridBlock_lowerOrder

    type :: gridBlock_t

        ! ---- index bookkeeping ----
        integer :: NGC = 0                                   ! guard cells padding each side
        integer, dimension(ndim) :: N = 0                    ! number of interior cells
        integer, dimension(ndim) :: minIdx = 0, maxIdx = 0   ! index of first/last guard cell
        integer, dimension(ndim) :: strtIdx = 0, stopIdx = 0 ! index of first/last interior cell

        ! ---- geometry ----
        real, dimension(ndim) :: domainBeg = 0.0, domainEnd = 0.0 ! bounds of this block
        real, dimension(ndim) :: dl = 0.0                         ! holds dx and dy
        real :: mindxdy = 0.0
        real :: maxdxdy = 0.0
        real, allocatable :: x(:), y(:)                           ! cell-center coordinates

        ! ---- pointwise face quadrature values of conservative variables ----
        ! Pointwise conservative variables reconstructed at face quadrature points
        ! (variable, quadrature_point, NSEW, xcoordinate, ycoordinate)
        real, allocatable :: faceVals(:,:,:,:,:)

        ! ---- face averaged fluxes ----
        ! face averaged fluxes found by taking weighted sum of the
        ! fluxes from each riemann problem that was solved at a quadrature point
        ! (variable, dir, xcoordinate, ycoordinate)
        real, allocatable :: flux(:,:,:,:)

        ! ---- MOOD scheme cascade ----
        ! which reconstruction scheme to use at each face of each cell
        ! 3->7th order GP-R3
        ! 2->5th order GP-R2
        ! 1->3rd order GP-R1
        ! 0->1st order Godunov (FOG)
        ! (faceIdx, i, j)
        integer, allocatable :: scheme(:,:,:)

        ! whether a cell's face values, fluxes, and RK substage update need to be
        ! recomputed on the next iteration of the MOOD loop inside Forward Euler.
        ! Set for a cell whose order was lowered and its four neighbors since they share
        ! faces with it
        ! (i, j)
        logical, allocatable :: recompute(:,:)

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

        ! face reconstructions 
        allocate(blk%faceVals(nConsVars, nquad, nfaces, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! fluxes 
        ! note that
        ! the flux at the (i-1/2,j) interface
        ! is given by flux(:,xdir,i,j)
        ! and 
        ! the flux at the (i+1/2,j) interface
        ! is given by flux(:,xdir,i+1,j)
        ! and 
        ! the flux at the (i,j-1/2) interface
        ! is given by flux(:,xdir,i,j-1)
        ! and 
        ! the flux at the (i,j+1/2) interface
        ! is given by flux(:,xdir,i,j+1)
        allocate(blk%flux(nConsVars, ndim, &
                          blk%minIdx(xdir):blk%maxIdx(xdir), &
                          blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! MOOD scheme/GP radius used at each face of each cell
        allocate(blk%scheme(nfaces, &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! MOOD recompute flags for each cell
        allocate(blk%recompute( &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))

        !!! zero all of these out !!!
        blk%x = 0.0
        blk%y = 0.0
        blk%faceVals = 0.0
        blk%flux = 0.0
        blk%scheme = 0
        blk%recompute = .false.

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
        deallocate(blk%faceVals)
        deallocate(blk%flux)
        deallocate(blk%scheme)
        deallocate(blk%recompute)

    end subroutine gridBlock_dealloc

    pure subroutine gridBlock_lowerOrder(blk, i, j)
        ! purpose:      Lower the order of the reconstruction used in cell
        !               (i,j) by one, down to first order, and flag the
        !               cell and its four neighbors to be recomputed in
        !               the MOOD loop. The neighbors are flagged since the
        !               faces they share with (i,j) are reconstructed at
        !               the lower order of the two cells (EDP1)
        !
        ! Inputs:       - blk (gridBlock_t) the block containing the cell
        !               - i, j (integer) the indices of the cell
        !
        ! Outputs:      - blk (gridBlock_t) with blk%scheme and
        !                 blk%recompute updated
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        integer, intent(in) :: i, j

        blk%scheme(:,i,j) = MAX(blk%scheme(:,i,j) - 1, 0)

        blk%recompute(i,   j  ) = .true.
        blk%recompute(i+1, j  ) = .true.
        blk%recompute(i-1, j  ) = .true.
        blk%recompute(i,   j+1) = .true.
        blk%recompute(i,   j-1) = .true.

    end subroutine gridBlock_lowerOrder

end module gridBlock
