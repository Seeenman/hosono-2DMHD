module detect

    use, intrinsic :: ieee_arithmetic, only: IEEE_IS_NAN, IEEE_IS_FINITE

    use definitions, only:

    implicit none

    private

contains

    function detect_PAD(dens, pres) result(admissible)
        ! physical admisibility condition
        implicit none
        real, intent(in) :: dens, pres
        logical :: admissible

        admissible=.false.

        if ((dens>0) and (pres>0)) admissible=.true.
            
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
            (.not. IEEE_IS_FINITE(pres)) then
            admissible=.false.
        end if
            
    end function detect_CAD

    function detect_DMP(var, nbhd, dx, dy) result(admissible)
        ! Discrete maximum principle
        implicit none
        implicit none
        real, intent(in) :: var, dx, dy
        real, intent(in) :: nbhd(5)
        logical :: admissible
        ! local variables
        real :: nbhd_max, nbhd_min

        admissible = .false.

        nbhd_max = MAXVAL(nbhd)
        nbhd_min = MINVAL(nbhd)

        ! check for micro-oscillations (plateau detection)
        if (nbhd_max-nbhd_min > min(dx,dy)) then
            ! no micro-oscillations now check DMP
            if ((MINVAL(nbhd)<=var) .and. (MAXVAL(nbhd)>=var)) then
                ! DMP is satisfied solution is admissible
                admissible = .true.
            else
                ! DMP is not satisfied, make solution admissible if u2 check passes
                
            end if
        end if

        ! u2 detection
        
        
    end function detect_DMP
    
end module detect
