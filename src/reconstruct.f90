module reconstruct

    use definitions, only: ndim, xdir, ydir, nConsVars
    use gridBlock, only: gridBlock_t
    use GP, only: GP_stencIdxs, GP_nStenc, GP_nPred, GP_predictionVectors, &
        GP_nQuadrature, GP_maxRadius

    implicit none

    private

    public :: reconstruct_faceValsWithGP

contains
    
    pure subroutine reconstruct_faceValsWithGP(blk, radius)
        ! purpose:   Reconstruct high order pointwise values at cell faces to be
        !            used to compute high order fluxes at cell faces. The order
        !            of the reconstruction is 2*radius+1.
        !            if a radius of 0 is given then this reverts to FOG.
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - The GP radius to use for this reconstruction which
        !              determines order of accuracy of the reconstruction as
        !              well as the number of quadrature points
        !            
        ! Outputs:   - blk%upperFace (real array) the Riemann state on the upper face
        !              of a cell in each direction (x and y) of each cell
        !            - blk%lowerFace (real array) the Riemann state upper face
        !              of a cell in each direction (x and y) of each cell
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        integer, intent(in) :: radius
        ! local variables
        integer :: ns 
        integer :: np 
        integer :: nq 
        real, allocatable :: stencil_data(:,:)
        real, allocatable :: faceVals(:,:)
        integer :: i,j
        integer :: ip,jp,k

        if (radius==0) then
            blk%upperFace(:, ydir, 1, :,:) = blk%U(:,:,:)
            blk%lowerFace(:, ydir, 1, :,:) = blk%U(:,:,:)
            blk%upperFace(:, xdir, 1, :,:) = blk%U(:,:,:)
            blk%lowerFace(:, xdir, 1, :,:) = blk%U(:,:,:)
        else if (radius > GP_maxRadius) then
            error stop ": GP radius larger than 3 not currently supported"
        else
            allocate(stencil_data(GP_nStenc(radius), nConsVars))
            allocate(faceVals(GP_nPred(radius), nConsVars))
            ns = GP_nStenc(radius)
            np = GP_nPred(radius)
            nq = GP_nQuadrature(radius)
            
            do j = blk%strtIdx(ydir)-1, blk%stopIdx(ydir)+1
                do i = blk%strtIdx(xdir)-1, blk%stopIdx(xdir)+1
                    do k=1, ns
                        ip = GP_stencIdxs(xdir, k, radius) + i
                        jp = GP_stencIdxs(ydir, k, radius) + j
                        stencil_data(k,:) = blk%U(:,ip,jp)
                    end do
                    faceVals = MATMUL(TRANSPOSE(GP_predictionVectors(1:ns, 1:np, radius)), stencil_data)
                    blk%upperFace(:, ydir, 1:nq, i,j) = TRANSPOSE(faceVals(0*nq+1:1*nq, :))
                    blk%lowerFace(:, ydir, 1:nq, i,j) = TRANSPOSE(faceVals(1*nq+1:2*nq, :))
                    blk%upperFace(:, xdir, 1:nq, i,j) = TRANSPOSE(faceVals(2*nq+1:3*nq, :))
                    blk%lowerFace(:, xdir, 1:nq, i,j) = TRANSPOSE(faceVals(3*nq+1:4*nq, :))
                end do
            end do
            deallocate(stencil_data,faceVals)
        end if


    end subroutine reconstruct_faceValsWithGP

end module reconstruct
