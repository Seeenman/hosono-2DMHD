module timeStep

    use definitions, only: xdir, ydir, nConsVars
    use gridBlock, only: gridBlock_t
    use gridState, only: gridState_t
    use convert, only: convert_cons2prim
    use boundaryConditions, only: boundaryConditions_applyBlk
    use reconstruct, only: reconstruct_faceValsWithGP
    use getFluxes, only: getFluxes_

    implicit none

    private

    public timeStep_advanceSolution

contains

    subroutine timeStep_advanceSolution(dt, blk, state)
        implicit none
        real, intent(in) :: dt
        type(gridBlock_t), intent(in out) :: blk
        type(gridState_t), intent(in out) :: state
        ! local variables
        integer :: i, j

        state%U = timeStep_forwardEuler(dt, blk, state%U)
        
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                state%V(:,i,j) = convert_cons2prim(state%U(:,i,j))
            end do
        end do

        call boundaryConditions_applyBlk(blk, state%U)
        call boundaryConditions_applyBlk(blk, state%V)
        
    end subroutine timeStep_advanceSolution

    function timeStep_forwardEuler(dt, blk, U) result(Unew)
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
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        real :: Unew(nConsVars, &
                     blk%minIdx(xdir):blk%maxIdx(xdir), &
                     blk%minIdx(ydir):blk%maxIdx(ydir))
        ! local variables
        real :: dx, dy
        integer i, j

        dx = blk%dl(xdir)
        dy = blk%dl(ydir)

        ! reconstruct face values at all quadrature points at the correct
        ! order based on blk%scheme. blk%scheme holds the radius of the
        ! GP method to be used for each face of each cell. 
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

    end function timeStep_forwardEuler

end module timeStep
