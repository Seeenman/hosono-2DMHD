module GP_volAvgToVolAvgSqExpKernel

    ! volume averaged to volume averaged squared exponential Gaussian Process kerenel function.

    use definitions, only: ndim, qp, qp_pi

    implicit none

contains

    function GP_volAvgToVolAvgSqExpKernel_(p1, p2, dl, ell) result(kernel)
        ! function:     GP_volAvgToVolAvgSqExpKernel_
        ! purpose:      Returns the covariance between volume averaged data located at p1=(x1,y1) and 
        !               volume averaged data located at p2=(x2,y2) where 
        !               we assume that pointwise data follows a Gaussian Process with a covariance function
        !               defined by the squared exponential (SE) kernel k(a,b). In other words, this
        !               function returns the integral of k(a,b) over both coordinates (a and b) where
        !               the integration in the first coordinate (the a-coordinate) of k is taken over the volume
        !               [x1-dx, x1+dx]x[y1-dy,y1+dy] and the integration in the second coordinate (the b-coordinate)
        !               of k is taken over the volume [x2-dx, x2+dx]x[y2-dy,y2+dy].
        !               
        !               See Reyes et al. 2018 "A New Class of High-Order Methods for Fluid Dynamics
        !               Simulations Using Gaussian Process Modeling: One-Dimensional Case" equation
        !               equation 22 or Bourgeois and Lee 2022 "GP-MOOD: A positivity-preserving high-order
        !               finite volume method for hyperbolic conservation laws" equations 9 and 12.
        ! 
        ! Inputs:       - p1(ndim) (quadruple precision real) spatial coordinate of the first data point, i.e.,
        !                 p1=(x1,y1)
        !               - p2(ndim) (quadruple precision real) spatial coordinate of the second data point, i.e.,
        !                 p2=(x2,y2)
        !               - dl(ndim) (quadruple precision real) grid spacing in x and y directions, i.e.,
        !                 dl(1) = dx and dl(2) = dy
        !               - ell (quadruple precision real) 
        !               
        ! Outputs:      - kernel (quadruple precision real) the value of the integrated SE kernel
        !                 evaluated at (p1, p2).
        ! ------------------------------------------------------------

        implicit none
        
        real(qp), intent(in), dimension(ndim) :: p1, p2, dl
        real(qp), intent(in) :: ell

        ! local variables
        real(qp) :: delta
        real(qp) :: kernel,
        real(qp) :: r1, r2, r3, r4, r5, r6
        real(qp) :: ell_over_dl
        integer :: i
        
        kernel = 1._qp
        do i=1,ndim
            delta = (p1(i)-p2(i))/dl(i) ! equation 13 in Bourgeois and Lee
            ell_over_dl = ell/dl(i)
            r1 = (delta+1._qp)/(SQRT(2._qp)*ell_over_dl)
            r2 = (delta-1._qp)/(SQRT(2._qp)*ell_over_dl)
            r3 = - (delta+1._qp)**2 / (2._qp*ell_over_dl**2)
            r4 = - (delta-1._qp)**2 / (2._qp*ell_over_dl**2)
            r5 = delta / (SQRT(2._qp)*ell_over_dl)
            r6 = delta**2 / (2._qp*ell_over_dl**2)

            kernel = kernel*&
                SQRT(qp_pi)*(ell_over_dl)**2 * &
                (&
                r1*ERF(r1) + r2*ERF(r2) &
                + 1._qp/SQRT(qp_pi)*(EXP(r3) + EXP(r4)) &
                - 2._qp*(r5*ERF(rf) + 1/SQRT(qp_pi)*EXP(r6)) &
                )
        end do


        return
    end function GP_volAvgToVolAvgSqExpKernel_

end module GP_volAvgToVolAvgSqExpKernel
