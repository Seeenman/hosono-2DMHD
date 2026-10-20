module getFluxes

    use definitions, only: nConsVars, xdir, ydir, north, south, east, west
    use gridBlock, only: gridBlock_t
    use riemannsolver, only: riemannsolver_getFaceFlux
    use GP, only: GP_nQuadratureMax
    implicit none

    private

    public :: getFluxes_

contains

    ! pure subroutine getFluxes_(blk)
    subroutine getFluxes_(blk)
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        ! Riemann states at every quadrature point on a face
        ! (variable, quadrature_point)
        real :: uL(nConsVars, GP_nQuadratureMax), uR(nConsVars, GP_nQuadratureMax)
        integer :: i, j
        integer :: rr


        !$omp parallel do collapse(2) default(none) shared(blk) private(rr,uL,uR)
        do j = blk%strtIdx(ydir)-1, blk%stopIdx(ydir)+2
            do i = blk%strtIdx(xdir)-1, blk%stopIdx(xdir)+2
                ! only recompute fluxes at faces of cells flagged by the
                ! MOOD loop. The face values on both sides of an unflagged
                ! face are unchanged so its flux is too.
                ! As this is currently writtin it will sometimes recompute
                ! fluxes that don't need to be recomputed but
                ! this should only be a marginal slow down
                ! since so few cells are flagged for recompute
                if (blk%recompute(i,j) .or. blk%recompute(i-1,j)) then
                    ! |             |             |            |
                    ! |     i-1   uL|uR    i      |     i+1    |
                    ! |             |             |            |
                    ! right riemann state at i-1/2 interface
                    uR = blk%faceVals(:, :, west, i, j)
                    ! left riemann state at i-1/2 interface
                    uL = blk%faceVals(:, :, east, i-1, j)

                    ! the GP radius used for this face
                    rr = MIN(blk%scheme(west,i,j), blk%scheme(east,i-1,j))

                    ! the flux at the i-1/2 interface
                    blk%fluxRaw(:, xdir, i, j) = riemannsolver_getFaceFlux(uL, uR, xdir, rr)
                    blk%flux(:, xdir, i, j) = blk%fluxRaw(:, xdir, i, j)
                end if


                ! ----------------
                ! 
                ! 
                !       j+1
                ! 
                ! 
                ! ----------------
                ! 
                ! 
                !       j
                ! 
                !       uR
                ! ----------------
                !       uL
                ! 
                !       j-1
                ! 
                ! 
                ! ----------------
                if (blk%recompute(i,j) .or. blk%recompute(i,j-1)) then
                    ! right riemann state at j-1/2 interface
                    uR = blk%faceVals(:, :, south, i, j)
                    ! left riemann state at i-1/2 interface
                    uL = blk%faceVals(:, :, north, i, j-1)

                    ! the GP radius used for this face
                    rr = MIN(blk%scheme(south,i,j), blk%scheme(north,i,j-1))

                    ! the flux at the j-1/2 interface
                    blk%fluxRaw(:, ydir, i, j) = riemannsolver_getFaceFlux(uL, uR, ydir, rr)
                    blk%flux(:, ydir, i, j) = blk%fluxRaw(:, ydir, i, j)
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine getFluxes_

end module getFluxes
