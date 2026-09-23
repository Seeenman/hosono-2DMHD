module reconstruct

    use definitions, only: ndim, xdir, ydir, nConsVars
    use gridBlock, only: gridBlock_t
    use GP, only: GP_stencIdxs, GP_nStenc, GP_nPred, GP_predictionVectors, &
        GP_nQuadrature

    implicit none

    private

contains
    
    pure subroutine reconstruct_faceValsWithGP(blk, radius)
        ! purpose:   Reconstruct high order pointwise values at cell faces to be
        !            used to compute high order fluxes at cell faces. The order
        !            of the reconstruction is 2*radius+1
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - The GP radius to use for this reconstruction which
        !              determines order of accuracy of the reconstruction as
        !              well as the number of quadrature points
        !            
        ! Outputs:   - blk%upperFace (real array) the Riemann States on the 
        !            - blk%lowerFace (real array) 
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        integer, intent(in) :: radius
        ! local variables
        integer :: ns = GP_nStenc(radius)
        integer :: np = GP_nPred(radius)
        integer :: nq = GP_nQuadrature(radius)
        real :: stencil_data(ns, nConsVars)
        real :: facVals(np, nConsVars)
        integer :: i,j
        integer :: ip,jp,k

        do j = blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i = blk%strtIdx(xdir), blk%stopIdx(xdir)
                do k=1, ns
                    ip = GP_stencIdxs(xdir, k, radius) + i
                    jp = GP_stencIdxs(ydir, k, radius) + j
                    stencil_data(k,:) = blk%U(:,ip,jp)
                end do
                faceVals = MATMUL(TRANSPOSE(GP_predictionVectors(1:ns, 1:np, radius)), stencil_data)
                blk%upperFace(:, ydir, 0*nq+1:1*nq, i,j) = facVals(0*nq+1:1*nq, :)
                blk%lowerFace(:, ydir, 1*nq+1:2*nq, i,j) = facVals(1*nq+1:2*nq, :)
                blk%upperFace(:, xdir, 2*nq+1:3*nq, i,j) = facVals(3*nq+1:3*nq, :)
                blk%lowerFace(:, xdir, 3*nq+1:4*nq, i,j) = facVals(3*nq+1:4*nq, :)
            end do
        end do


    end subroutine reconstruct_faceValsWithGP

end module reconstruct
