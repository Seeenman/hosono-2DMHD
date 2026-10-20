module FLICDdivB 
    use definitions, only: xdir, ydir, magx_var, magy_var
    use gridBlock, only: gridBlock_t

    implicit none

    private

    public :: FLICDdivB_fluxCorrections

contains

    subroutine FLICDdivB_fluxCorrections(blk)
        ! Flux-interpolated central difference scheme to 
        ! maintain div(B)=0 from Tóth section 4.5.

        implicit none
        type(gridBlock_t), intent(in out) :: blk
        ! local variables
        integer :: i,j

        ! raw flux holds the fluxes that the Riemann solver
        ! returns without any div(B)=0 control. rawFlux is used here
        ! to build the corrected components of the magnetic field
        ! by using the cell centered constrained transport scheme from
        ! Tóth's 2000 paper "the div(B)=0 constraint in Shock-Capturing MHD codes".
        ! This is necessary because otherwise this routine
        ! would overwrite blk%flux multiple times during the MOOD loop
        ! and we wouldn't have a consistent starting point to recompute the
        ! corrected values of the magnetic field from on each iteration
        ! of the mood loop.

        !$omp parallel do collapse(2) private(i,j)
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir) + 1
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir) + 1
                if (blk%recomputeFLICD(i,j)) then
                    ! correction for the magy_var component of the xdir flux
                    blk%flux(magy_var,xdir,i,j) = &
                        (2 * blk%fluxRaw(magy_var,xdir,i  ,j  ) &
                        + blk%fluxRaw(magy_var,xdir,i+1,j  ) & 
                        + blk%fluxRaw(magy_var,xdir,i-1,j  ) & 
                        - blk%fluxRaw(magx_var,ydir,i  ,j+1) & 
                        - blk%fluxRaw(magx_var,ydir,i  ,j  ) & 
                        - blk%fluxRaw(magx_var,ydir,i-1,j  ) & 
                        - blk%fluxRaw(magx_var,ydir,i-1,j+1) )/8.0
                    ! correction for the magx_var component of the ydir flux
                    blk%flux(magx_var,ydir,i,j) = &
                        (2 * blk%fluxRaw(magx_var,ydir,i  ,j  ) &
                        + blk%fluxRaw(magx_var,ydir,i  ,j+1) & 
                        + blk%fluxRaw(magx_var,ydir,i  ,j-1) & 
                        - blk%fluxRaw(magy_var,xdir,i+1,j  ) & 
                        - blk%fluxRaw(magy_var,xdir,i  ,j  ) & 
                        - blk%fluxRaw(magy_var,xdir,i  ,j-1) & 
                        - blk%fluxRaw(magy_var,xdir,i+1,j-1) )/8.0
                end if
            end do
        end do
        !$omp end parallel do

        ! The normal component of B must carry no flux of its own. In the 2D FLICD
        ! scheme B_x is advanced solely by the y derivative of the z-component of the 
        ! electric field (Ez) and B_y is advanced solely by the x derivative of Ez.
        ! So F^x(B_x) and F^y(B_y) have to vanish for the divergence terms to cancel.
        ! The analytic flux is already zero, but some Riemann solvers may (e.g. HLLC)
        ! add a dissipation that is not part of the cancellation. 
        ! So we zero out those components here.
        blk%flux(magx_var, xdir, :, :) = 0.0
        blk%flux(magy_var, ydir, :, :) = 0.0
        
    end subroutine FLICDdivB_fluxCorrections
        
end module FLICDdivB 
