module forwardEuler

    use definitions, only: ndim, xdir, ydir, nConsVars, north, south, east, west
    use gridBlock, only: gridBlock_t
    use reconstruct, only: reconstruct_faceValsWithGP
    use getFluxes, only: getFluxes_

    implicit none

    private

    public :: forwardEuler_

contains

    subroutine forwardEuler_(blk, U)
        ! purpose:   
        !            
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - U (real array) the conservative variables on blk
        !              to advance, e.g., grid_state%U or the U of an RK
        !              substage
        !            
        ! Outputs:   
        !            
        !            
        !            
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        ! local variables
        real :: dx, dy
        real, dimension(nConsVars) :: Fplus, Fminus, Gplus, Gminus
        integer i, j, f

        dx = blk%dl(xdir)
        dy = blk%dl(ydir)

        ! reconstruct face values at all quadrature points at the correct
        ! order based on blk%scheme. blk%scheme holds the radius of the
        ! GP method to be used for each face of each cell. 
        call reconstruct_faceValsWithGP(blk, U)

        ! fill in fluxes
        call getFluxes_(blk)

    end subroutine forwardEuler_

end module forwardEuler
