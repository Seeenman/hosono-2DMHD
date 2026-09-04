module cfl


    use definitions
    use simulation, only: sim_cfl, sim_forceHydro
    use eigen, only: eigen_c_fs
    use gridBlock, only: gridBlock_t

    implicit none

contains

    function cfl_computedt(blk) result(dt)
        ! Based on equation 16.38 in Toro. Can underestimate max wave speed

        type(gridBlock_t), intent(in) :: blk
        real :: dt
        ! local variables
        real :: xMaxSpeed, xLambda, &
                yMaxSpeed, yLambda, &
                soundSpeed, BB, B_N_min
        integer :: i, j

        xMaxSpeed = 0.0
        yMaxSpeed = 0.0
        !$omp parallel do collapse(2) private(i,j,soundSpeed,xLambda,yLambda,BB,B_N_min) reduction(max:xMaxSpeed,yMaxSpeed)
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                if (sim_forceHydro) then
                    soundSpeed = sqrt(blk%V(gamm_var,i,j)*blk%V(pres_var,i,j)/blk%V(dens_var,i,j))
                else
                    ! For MHD use fast magneto-acoustic speed instead of sound speed
                    BB = blk%V(magx_var,i,j)**2 + blk%V(magy_var,i,j)**2 + blk%V(magz_var,i,j)**2
                    B_N_min = MINVAL(ABS(blk%V(magx_var:magz_var,i,j))) ! is this line correct? TODO
                    soundSpeed = eigen_c_fs(blk%V(gamm_var,i,j), &
                        blk%V(dens_var,i,j), &
                        blk%V(pres_var,i,j), &
                        BB, B_N_min, +1.0)
                end if
                xLambda = abs(blk%V(velx_var,i,j)) + soundSpeed
                yLambda = abs(blk%V(vely_var,i,j)) + soundSpeed
                xMaxSpeed = max(xMaxSpeed, xLambda) 
                yMaxSpeed = max(yMaxSpeed, yLambda) 
            end do
        end do
        !$omp end parallel do

        dt = sim_cfl * min(blk%dl(XDIR)/xMaxSpeed, blk%dl(YDIR)/yMaxSpeed)

    end function cfl_computedt

end module cfl
