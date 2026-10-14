module detect

    use, intrinsic :: ieee_arithmetic, only: IEEE_IS_NAN, IEEE_IS_FINITE

    use definitions, only: nConsVars, nPrimVars, dens_var, pres_var, velx_var, vely_var, pres_var, gamm_var, &
        xdir, ydir, ndim
    use gridBlock, only: gridBlock_t
    use convert, only: convert_cons2prim

    implicit none

    private
    public :: detect_ALL
    ! public :: detect_CAD

contains

    subroutine detect_ALL(blk, U, Unew, exit_mood_loop)
        ! purpose:   TODO
        !            
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - U (real array) the conservative variables on blk
        !              from the last RK substage
        !            - U (real array) the candidate conservative variables on blk
        !              from the last RK substage
        !            
        ! Outputs:   TODO
        !            
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        real, intent(in out) :: Unew(nConsVars, &
                                     blk%minIdx(xdir):blk%maxIdx(xdir), &
                                     blk%minIdx(ydir):blk%maxIdx(ydir))
        logical, intent(out) :: exit_mood_loop
        ! local variables
        real :: dx, dy
        real :: div_vel, mach, pres_gradient
        real :: densNew, dens1, dens2, dens3, dens4, dens5, min_dens, max_dens
        real :: dens2ndDerivs(5,2), min_dens2ndDerivs(2), max_dens2ndDerivs(2)
        real :: sigma_v, sigma_M, sigma_p
        real, dimension(nPrimVars) :: V, Vi_plus, Vi_minus, Vj_plus, Vj_minus
        integer :: i, j
        integer :: ii, jj, counter, dir
        logical :: admissible, shock1, shock2, shock3
        logical :: sx1, sx2, sx3, sy1, sy2, sy3, smooth1, smooth2

        exit_mood_loop = .true.

        dx = blk%dl(xdir)
        dy = blk%dl(ydir)

        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)

                ! -------------------------
                ! check for NANs and INFs
                ! -------------------------
                admissible = detect_CAD(Unew(dens_var,i,j), Unew(pres_var, i, j))
                if (.not. admissible) then
                    blk%scheme(:,i,j) = blk%scheme(:,i,j) - 1
                    exit_mood_loop = .false.
                    CYCLE ! go to next index
                end if

                ! Check for positivity
                admissible = detect_PAD(Unew(dens_var,i,j), Unew(pres_var, i, j))
                if (.not. admissible) then
                    blk%scheme(:,i,j) = blk%scheme(:,i,j) - 1
                    exit_mood_loop = .false.
                    CYCLE ! go to next index
                end if

                ! TODO move this check into the mood loop because it only needs
                ! to be computed once per RK substage
                ! -------------------------
                ! check for shocks
                ! -------------------------
                shock1 = .false.
                shock2 = .false.
                shock3 = .false.
                V = convert_cons2prim(U(:,i,j))
                Vi_plus = convert_cons2prim(U(:,i+1,j))
                Vi_minus = convert_cons2prim(U(:,i-1,j))
                Vj_plus = convert_cons2prim(U(:,i,j+1))
                Vj_minus = convert_cons2prim(U(:,i,j-1))

                ! shock check 1: local compressibility
                div_vel = (Vi_plus(velx_var) - Vi_minus(velx_var))/(2*dx)
                div_vel = div_vel + (Vj_plus(vely_var) - Vj_minus(vely_var))/(2*dy)
                sigma_v = 1.0 ! TODO make sigma_v a paramter that is read in from the par file
                if (div_vel < - sigma_v*blk%maxdxdy**2) shock1 = .true. ! TODO make sigma_v a paramter that is read in from the par file

                ! shock check 2: mach number
                mach = V(velx_var)**2 + V(vely_var)**2
                mach = V(dens_var)*mach/(V(gamm_var)*V(pres_var))
                mach = SQRT(mach)
                sigma_M = 0.2 ! TODO make sigma_M a paramter that is read in from the par file
                if (mach > sigma_M) shock2 = .true.

                ! shock check 3: normalized local pressure gradient
                pres_gradient = ABS(Vi_plus(pres_var) - Vi_minus(pres_var))/(MIN(Vi_plus(pres_var), Vi_minus(pres_var)))
                pres_gradient = pres_gradient + ABS(Vj_plus(pres_var) - Vj_minus(pres_var))/(MIN(Vj_plus(pres_var), Vj_minus(pres_var)))
                sigma_p = 10.0 ! TODO make sigma_p a paramter that is read in from the par file
                if (pres_gradient > sigma_p*blk%maxdxdy) shock3 = .true.

                ! if any of the shock tests return false
                if (.not. (shock1 .and. shock2 .and. shock3)) then
                    CYCLE ! go to next index
                end if

                ! --------------------------------------------
                ! plateau, discrete maximum principle, and u2 
                ! --------------------------------------------
                densNew = Unew(dens_var, i,j)
                dens1 = U(dens_var, i,j)
                dens2 = U(dens_var, i+1,j)
                dens3 = U(dens_var, i-1,j)
                dens4 = U(dens_var, i,j+1)
                dens5 = U(dens_var, i,j-1)
                min_dens = MIN(dens1, dens2, dens3, dens4, dens5)
                max_dens = MAX(dens1, dens2, dens3, dens4, dens5)

                ! check for micro-oscillations (plateau detection)
                if (max_dens - min_dens < blk%mindxdy**3) then
                    ! micro-oscillations are okay
                    CYCLE ! go to next index
                end if

                ! Check DMP if no micro-oscillations
                if ((min_dens <= densNew) .and. (densNew <= max_dens)) then
                    ! DMP met
                    CYCLE ! go to next index
                end if

                ! TODO this can also be moved into the MOOD loop and only computed once
                ! per RK substage and checked before plateau and DMP. But it
                ! probably isn't worth the effort since this check
                ! almost never runs in practice I think
                ! u2 detection  think
                counter = 0
                do ii=-1,1; do jj = -1, 1
                    if (abs(ii) + abs(jj) <= 1) then
                        counter = counter + 1
                        dens2ndDerivs(counter,xdir) = (U(dens_var, i+ii+1, j+jj) - 2*U(dens_var, i+ii, j+jj) + U(dens_var, i+ii-1, j+jj))/dx**2
                        dens2ndDerivs(counter,ydir) = (U(dens_var, i+ii, j+jj+1) - 2*U(dens_var, i+ii, j+jj) + U(dens_var, i+ii, j+jj-1))/dy**2
                    end if
                end do; end do

                min_dens2ndDerivs(xdir) = MINVAL(dens2ndDerivs(:,xdir))
                min_dens2ndDerivs(ydir) = MINVAL(dens2ndDerivs(:,ydir))
                max_dens2ndDerivs(xdir) = MAXVAL(dens2ndDerivs(:,xdir))
                max_dens2ndDerivs(ydir) = MAXVAL(dens2ndDerivs(:,ydir))

                sx1 = (min_dens2ndDerivs(xdir)*max_dens2ndDerivs(xdir) > -blk%mindxdy)
                sy1 = (min_dens2ndDerivs(ydir)*max_dens2ndDerivs(ydir) > -blk%mindxdy)
                sx2 = (MAX(ABS(min_dens2ndDerivs(xdir)), ABS(max_dens2ndDerivs(xdir))) < blk%mindxdy)
                sy2 = (MAX(ABS(min_dens2ndDerivs(ydir)), ABS(max_dens2ndDerivs(ydir))) < blk%mindxdy) 
                sx3 = ((ABS(min_dens2ndDerivs(xdir))/ABS(max_dens2ndDerivs(xdir))) > 0.5)
                sy3 = ((ABS(min_dens2ndDerivs(ydir))/ABS(max_dens2ndDerivs(ydir))) > 0.5)

                smooth1 = (sx1 .and. sx2 .and. sy1 .and. sy2)
                smooth2 = (sx1 .and. sx3 .and. sy1 .and. sy3)

                if (.not. (smooth1 .or. smooth2)) then
                    blk%scheme(:,i,j) = blk%scheme(:,i,j) - 1
                    exit_mood_loop = .false.
                end if
                
            end do
        end do
        
    end subroutine detect_ALL

    function detect_PAD(dens, pres) result(admissible)
        ! physical admisibility condition
        implicit none
        real, intent(in) :: dens, pres
        logical :: admissible

        admissible=.false.

        if ((dens>0.0) .and. (pres>0.0)) admissible=.true.
            
    end function detect_PAD

    function detect_CAD(dens, pres) result(admissible)
        ! Computer science admissibility detection
        implicit none
        real, intent(in) :: dens, pres
        logical :: admissible

        admissible=.true.

        if ((IEEE_IS_NAN(dens)) .or. &
            (IEEE_IS_NAN(pres)) .or. &
            (.not. IEEE_IS_FINITE(dens)) .or. &
            (.not. IEEE_IS_FINITE(pres))) then
            admissible=.false.
        end if
            
    end function detect_CAD
    
end module detect
