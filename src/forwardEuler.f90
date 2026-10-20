module forwardEuler

    use definitions, only: ndim, xdir, ydir, nConsVars, north, south, east, west
    use gridBlock, only: gridBlock_t
    use reconstruct, only: reconstruct_faceValsWithGP
    use getFluxes, only: getFluxes_

    implicit none

    private

    public :: forwardEuler_

contains

    function forwardEuler_(dt, blk, U) result(Unew)
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
        real, intent(in) :: dt
        type(gridBlock_t), intent(in) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        real :: Unew(nConsVars, &
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
        ! GP method to be used in each cell. 
        call reconstruct_faceValsWithGP(blk, U)

        ! fill in fluxes. Don't need to pass in U because this reads
        ! directly from blk%faceVals which was updated by
        ! reconstruct_faceValsWithGP using whatever substage value of U
        ! we are on.
        call getFluxes_(blk)

        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                Unew(:,i,j) = U(:,i,j) &
                    - (dt/dx)*(blk%flux(:,xdir,i+1,j) - blk%flux(:,xdir,i,j)) &
                    - (dt/dy)*(blk%flux(:,ydir,i,j+1) - blk%flux(:,ydir,i,j))
            end do
        end do

    end function forwardEuler_

end module forwardEuler
