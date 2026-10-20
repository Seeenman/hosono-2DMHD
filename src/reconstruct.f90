module reconstruct

    use definitions, only: ndim, xdir, ydir, nConsVars, dens_var, north, south, east, west
    use gridBlock, only: gridBlock_t, gridBlock_lowerOrder
    use GP, only: GP_stencIdxs, GP_nStenc, GP_nPred, GP_predictionVectors, &
        GP_nQuadrature, GP_maxRadius, GP_nPredMax, GP_nQuadratureMax, GP_nStencMax
    use detect, only: detect_CAD, detect_PAD
    use convert, only: convert_cons2pres

    implicit none

    private

    public :: reconstruct_faceValsWithGP

contains

    subroutine reconstruct_faceValsWithGP(blk, U)
        ! purpose:   Reconstruct high order pointwise values at cell faces to be
        !            used to compute high order fluxes at cell faces. 
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - U (real array) the conservative variables (cell
        !              averages) on blk to reconstruct, e.g., grid_state%U
        !              or the U of an RK substage
        !            
        ! Outputs:   - blk%upperFace (real array) the Riemann state on the upper face
        !              of a cell in each direction (x and y) of each cell
        !            - blk%lowerFace (real array) the Riemann state upper face
        !              of a cell in each direction (x and y) of each cell
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        ! local variables
        integer :: ns 
        integer :: nq 
        real :: stencil_data(nConsVars,GP_nStencMax)
        integer :: rr
        integer :: i,j
        integer :: f ! face index (1=north, 2=south, 3=east, 4=west)
        integer :: ia, ja ! offset to the neighbor sharing face f
        integer :: ip,jp,k
        integer :: q
        real :: pres
        logical :: faceValsAreAdmissible
        logical, allocatable :: cellIsAdmissible(:,:), reconstructCell(:,:), reconstructNextPass(:,:)


        ! reconstruct on the interior plus two rings of guard cells. getFluxes_
        ! needs face values one cell beyond the interior to get the fluxes on
        ! the boundary faces of the interior, plus one more cell to get the
        ! fluxes one face further out (needed for cell centered constrained
        ! transport). Requires NGC >= GP_maxRadius+2 (see grid.f90).
        ! only reconstruct cells flagged by the MOOD loop. Their neighbors'
        ! faces are unchanged
        allocate(cellIsAdmissible(blk%minIdx(xdir):blk%maxIdx(xdir), &
                                  blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(reconstructCell(blk%minIdx(xdir):blk%maxIdx(xdir), &
                                 blk%minIdx(ydir):blk%maxIdx(ydir)))
        allocate(reconstructNextPass(blk%minIdx(xdir):blk%maxIdx(xdir), &
                                     blk%minIdx(ydir):blk%maxIdx(ydir)))
        reconstructCell = blk%recompute

        faceValsAreAdmissible = .false.

        ! reconstruct until all Riemann states are admissible. First order
        ! face values are the cell averages, so this loop always terminates.
        ! First order values are not checked for admissibility.
        do while (.not. faceValsAreAdmissible)

            cellIsAdmissible = .true.

            !$omp parallel do collapse(2) default(none) &
            !$omp shared(blk, U, GP_nStenc, GP_nQuadrature, GP_predictionVectors, &
            !$omp GP_stencIdxs, cellIsAdmissible, reconstructCell) &
            !$omp private(k,f,ia,ja,ip,jp,rr,q,ns,nq,pres,stencil_data)
            do j = blk%strtIdx(ydir)-2, blk%stopIdx(ydir)+2
                do i = blk%strtIdx(xdir)-2, blk%stopIdx(xdir)+2
                    if (.not. reconstructCell(i,j)) CYCLE

                    do f=1,4 ! loop over all four faces of the cell

                        if (f==north) then
                            ia=0; ja=1
                        else if (f==south) then
                            ia=0; ja=-1
                        else if (f==east) then
                            ia=1; ja=0
                        else if (f==west) then
                            ia=-1; ja=0
                        end if

                        ! get the GP radius that we are reconstructing this face with
                        ! according to the EDP1 method from the GP-MOOD paper
                        rr = MIN(blk%scheme(i,j), blk%scheme(i+ia,j+ja))

                        if (rr==0) then ! FOG
                            blk%faceVals(1:nConsVars, 1, f, i, j) = U(1:nConsVars,i,j)
                        else ! high order GP reconstruction
                            ! get nStenc, and nQuadrature corresponding to the order at which
                            ! we will be reconstructing at this cell face
                            ns = GP_nStenc(     rr)
                            nq = GP_nQuadrature(rr)
                        
                            do k=1, ns
                                ip = GP_stencIdxs(xdir, k, rr) + i
                                jp = GP_stencIdxs(ydir, k, rr) + j
                                stencil_data(1:nConsVars,k) = U(1:nConsVars,ip,jp)
                            end do

                            ! Get the pointwise values at face quadrature points.
                            ! Written as explicit loops rather than with MATMUL and TRANSPOSE
                            ! to avoid allocating array temporaries on the heap for every face
                            do q=1, nq
                                blk%faceVals(1:nConsVars, q, f, i, j) = 0.0
                                do k=1, ns
                                    blk%faceVals(1:nConsVars, q, f, i, j) = blk%faceVals(1:nConsVars, q, f, i, j) + &
                                        GP_predictionVectors(k, q, f, rr)*stencil_data(1:nConsVars, k)
                                end do
                            end do

                            ! CAD and PAD check on the Riemann states
                            do q=1, nq
                                pres = convert_cons2pres(blk%faceVals(1:nConsVars, q, f, i, j))
                                if (.not. (detect_CAD(blk%faceVals(dens_var, q, f, i, j), pres) .and. &
                                           detect_PAD(blk%faceVals(dens_var, q, f, i, j), pres))) then
                                    cellIsAdmissible(i,j) = .false.
                                end if
                            end do

                        end if ! matches: else ! high order GP reconstruction

                    end do ! matches: do f=1,4 ! loop over all four faces of the cell

                end do ! matches: do i = blk%strtIdx(xdir)-2, blk%stopIdx(xdir)+2
            end do ! matches: do j = blk%strtIdx(ydir)-2, blk%stopIdx(ydir)+2
            !$omp end parallel do

            ! This has to be a separate serial loop because we can't
            ! write to neighboring cells inside of a loop wrapped in an
            ! openmp directive
            faceValsAreAdmissible = .true.
            reconstructNextPass = .false.
            do j = blk%strtIdx(ydir)-2, blk%stopIdx(ydir)+2
                do i = blk%strtIdx(xdir)-2, blk%stopIdx(xdir)+2
                    ! lower the order of a cell with inadmissible Riemann states.
                    ! Since faces use the lower order of the two cells sharing
                    ! them (EDP1 from GP-MOOD paper), the neighbors' faces are lowered
                    ! as well so both Riemann states at every face stay at the same order.
                    ! gridBlock_lowerOrder lowers the order of the reconstruction at
                    ! all four faces of cell (i,j). gridBlock_lowerOrder also sets
                    ! blk%recompute to .true. for cell (i,j) (which is the cell with
                    ! inadmissible Riemann states in this case) as well as its 4 directly
                    ! adjacent neighbors. 
                    ! This is done so that their fluxes and cell average updates are
                    ! recomputed and they are checked by the MOOD loop.
                    if (.not. cellIsAdmissible(i,j)) then
                        call gridBlock_lowerOrder(blk, i, j)
                        reconstructNextPass(i-1:i+1, j) = .true.
                        reconstructNextPass(i, j-1:j+1) = .true.
                        faceValsAreAdmissible = .false.
                    end if

                end do
            end do

            reconstructCell = reconstructNextPass

        end do

    end subroutine reconstruct_faceValsWithGP
    
end module reconstruct
