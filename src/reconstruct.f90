module reconstruct

    use definitions, only: ndim, xdir, ydir, nConsVars, north, south, east, west
    use gridBlock, only: gridBlock_t
    use GP, only: GP_stencIdxs, GP_nStenc, GP_nPred, GP_predictionVectors, &
        GP_nQuadrature, GP_maxRadius, GP_nPredMax, GP_nQuadratureMax, GP_nStencMax

    implicit none

    private

    public :: reconstruct_faceValsWithGP

contains

    pure subroutine reconstruct_faceValsWithGP(blk)
        ! purpose:   Reconstruct high order pointwise values at cell faces to be
        !            used to compute high order fluxes at cell faces. 
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            
        ! Outputs:   - blk%upperFace (real array) the Riemann state on the upper face
        !              of a cell in each direction (x and y) of each cell
        !            - blk%lowerFace (real array) the Riemann state upper face
        !              of a cell in each direction (x and y) of each cell
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        ! local variables
        integer :: ns 
        integer :: nq 
        real :: stencil_data(GP_nStencMax,nConsVars)
        integer :: rr
        integer :: i,j
        integer :: f ! face index (1=north, 2=south, 3=east, 4=west)
        integer :: fa, ia, ja ! adjacent face index and coordinates
        integer :: ip,jp,k


        do j = blk%strtIdx(ydir)-1, blk%stopIdx(ydir)+1
            do i = blk%strtIdx(xdir)-1, blk%stopIdx(xdir)+1
                do f=1,4 ! loop over all four faces of the cell

                    if (f==north) then
                        fa = south; ia=0; ja=1
                    else if (f==south) then
                        fa = north; ia=0; ja=-1
                    else if (f==east) then
                        fa = west; ia=1; ja=0
                    else if (f==west) then
                        fa = east; ia=-1; ja=0
                    end if

                    ! get the GP radius that we are reconstructing this face with
                    rr = MIN(blk%scheme(f,i,j), blk%scheme(fa,i+ia,j+ja))

                    ! get nStenc, and nQuadrature corresponding to the order at which
                    ! we will be reconstructing at this cell face
                    ns = GP_nStenc(     rr)
                    nq = GP_nQuadrature(rr)

                    do k=1, ns
                        ip = GP_stencIdxs(xdir, k, rr) + i
                        jp = GP_stencIdxs(ydir, k, rr) + j
                        stencil_data(k,1:nConsVars) = blk%U(1:nConsVars,ip,jp)
                    end do

                    ! Get the pointwise values at face quadrature points
                    blk%faceVals(1:nConsVars, 1:nq, f, i, j) = TRANSPOSE(MATMUL(&
                        TRANSPOSE(GP_predictionVectors(1:ns, 1:nq, f, rr)),&
                        stencil_data(1:ns, 1:nConsVars)))

                end do
            end do
        end do

    end subroutine reconstruct_faceValsWithGP
    
end module reconstruct
