module FLICDdivB 
    use definitions, only: xdir, ydir, ndim, magx_var, magy_var
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
        real :: newFlux(magx_var:magy_var, ndim, &
                        blk%minIdx(xdir):blk%maxIdx(xdir), &
                        blk%minIdx(ydir):blk%maxIdx(ydir))

        !$omp parallel do collapse(2) private(i,j)
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir) + 1
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir) + 1
                ! correction for the magy_var component of the xdir flux
                newFlux(magy_var,xdir,i,j) = &
                    2 * blk%flux(magy_var,xdir,i  ,j  ) &
                      + blk%flux(magy_var,xdir,i+1,j  ) & 
                      + blk%flux(magy_var,xdir,i-1,j  ) & 
                      - blk%flux(magx_var,ydir,i  ,j+1) & 
                      - blk%flux(magx_var,ydir,i  ,j  ) & 
                      - blk%flux(magx_var,ydir,i-1,j  ) & 
                      - blk%flux(magx_var,ydir,i-1,j  )
                newFlux(magy_var, xdir, i,j) = newFlux(magy_var, xdir, i,j)/8.0
                ! correction for the magx_var component of the ydir flux
                newFlux(magx_var,ydir,i,j) = &
                    2 * blk%flux(magx_var,ydir,i  ,j  ) &
                      + blk%flux(magx_var,ydir,i  ,j+1) & 
                      + blk%flux(magx_var,ydir,i  ,j-1) & 
                      - blk%flux(magy_var,xdir,i+1,j+1) & 
                      - blk%flux(magy_var,xdir,i  ,j  ) & 
                      - blk%flux(magy_var,xdir,i  ,j-1) & 
                      - blk%flux(magy_var,xdir,i+1,j-1)
                newFlux(magx_var, ydir, i,j) = newFlux(magx_var, ydir, i,j)/8.0
            end do
        end do
        !$omp end parallel do

        blk%flux(magy_var, xdir, :, :) = newFlux(magy_var, xdir, :,:)
        blk%flux(magx_var, ydir, :, :) = newFlux(magx_var, ydir, :,:)

        ! The normal component of B must carry no flux of its own. In the 2D FLICD
        ! scheme B_x is advanced solely by the y derivative of the z-component of the 
        ! electric field (Ez) and B_y is advanced solely by the x derivative of Ez.
        ! So F^x(B_x) and F^y(B_y) have to vanish for the divergence terms to cancel.
        ! The analytic flux is already zero, but some Riemann solvers may (e.g. HLLC)
        ! add a dissipative piece that is not part of the cancellation. 
        ! So we zero out those components here.
        blk%flux(magx_var, xdir, :, :) = 0.0
        blk%flux(magy_var, ydir, :, :) = 0.0
        
    end subroutine FLICDdivB_fluxCorrections
        
end module FLICDdivB 
