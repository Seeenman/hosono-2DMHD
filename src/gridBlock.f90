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
        real, allocatable :: fluxRaw(:,:,:,:)

        ! ---- MOOD scheme cascade ----
        ! which reconstruction scheme to use in each cell
        ! 3->7th order GP-R3
        ! 2->5th order GP-R2
        ! 1->3rd order GP-R1
        ! 0->1st order Godunov (FOG)
        ! A face is reconstructed with the lower scheme of the two cells
        ! sharing it (EDP1 from the GP-MOOD paper)
        ! (i, j)
        integer, allocatable :: scheme(:,:)

        ! whether a cell's face values, fluxes, and RK substage update need to be
        ! recomputed on the next iteration of the MOOD loop inside Forward Euler.
        ! Set for a cell whose order was lowered and its four neighbors since they share
        ! faces with it
        ! (i, j)
        logical, allocatable :: recompute(:,:)

        ! whether a cell's corrected fluxes (corrected by FLICD) get
        ! recomputed on the next iteration of the MOOD loop inside of Forward Euler.
        ! Also means that the unsplit finite volume update in timeStep_forwardEuler
        ! gets recomputed.
        ! Set for a cell whose order was lowered by detect_ALL and any cell within
        ! a radius of 2 under the taxi cab metric since if cell (i,j) get's its
        ! order lowered then all four face fluxes of that cell change and those
        ! fluxes are used by cells up to a distance of 2 away from cell (i,j)
        ! when the fluxes get corrected for cell centered constrained transport
        ! (i, j)
        logical, allocatable :: recomputeFLICD(:,:)

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
        ! raw flux holds the fluxes that the Riemann solver
        ! returns without any div(B)=0 control. fluxRaw is used by FLICDdivB
        ! to build the corrected components of the magnetic field
        ! by using the cell centered constrained transport scheme.
        ! This is necessary because otherwise FLICDdivB_fluxCorrections
        ! would overwrite blk%flux multiple times during the MOOD loop
        ! and we wouldn't have a consistent starting point to recompute the
        ! corrected values of the magnetic field from on each iteration
        ! of the mood loop.
        allocate(blk%fluxRaw(nConsVars, ndim, &
                             blk%minIdx(xdir):blk%maxIdx(xdir), &
                             blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! MOOD scheme/GP radius used in each cell
        allocate(blk%scheme( &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))

        ! MOOD recompute flags for each cell
        allocate(blk%recompute( &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(blk%recomputeFLICD( &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)))

        !!! zero all of these out !!!
        blk%x = 0.0
        blk%y = 0.0
        blk%faceVals = 0.0
        blk%flux = 0.0
        blk%fluxRaw = 0.0
        blk%scheme = 0
        blk%recompute = .false.
        blk%recomputeFLICD = .false.

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
        deallocate(blk%fluxRaw)
        deallocate(blk%scheme)
        deallocate(blk%recompute)
        deallocate(blk%recomputeFLICD)

    end subroutine gridBlock_dealloc

    pure subroutine gridBlock_lowerOrder(blk, i, j)
        ! purpose:      Lower the order of the reconstruction used in cell
        !               (i,j) by one, going no lower than first order, and flag the
        !               cell and its four neighbors to be recomputed in
        !               the MOOD loop. The neighbors are flagged since the
        !               faces they share with (i,j) are reconstructed at
        !               the lower order of the two cells (EDP1 from GP-MOOD paper).
        !               Additionally, all cells within radius 2 of (i,j) under
        !               the taxi cab metric are also flagged since cell centered
        !               constrained transport (for div(B)=0) requires an even larger
        !               stencil
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
        ! local variables
        integer :: ii, jj

        blk%scheme(i,j) = MAX(blk%scheme(i,j) - 1, 0)

        blk%recompute(i,   j  ) = .true.
        blk%recompute(i+1, j  ) = .true.
        blk%recompute(i-1, j  ) = .true.
        blk%recompute(i,   j+1) = .true.
        blk%recompute(i,   j-1) = .true.

        ! every cell within 2 of (i,j) in the taxicab metric, clipped to the
        ! array since reconstruct lowers cells as far out as strtIdx-2, which
        ! is only 1 cell from the edge of the array when GP_maxRadius=1
        do jj = MAX(j-2, blk%minIdx(ydir)), MIN(j+2, blk%maxIdx(ydir))
            do ii = MAX(i-2, blk%minIdx(xdir)), MIN(i+2, blk%maxIdx(xdir))
                if (ABS(ii-i) + ABS(jj-j) <= 2) blk%recomputeFLICD(ii,jj) = .true.
            end do
        end do

    end subroutine gridBlock_lowerOrder

end module gridBlock
