module GP_volAvgToPointSqExpKernel

    ! volume averaged to pointwise squared exponential Gaussian Process kerenel function.

    use definitions, only: ndim, qp, qp_pi

    implicit none

contains

    function GP_volAvgToPointSqExpKernel_(p1, p2, dl, ell) result(kernel)
        ! function:     GP_volAvgToPointSqExpKernel_
        ! purpose:      Returns the covariance between volume averaged data located at p1=(x1,y1) and 
        !               pointwise data located at p2=(x2,y2) where 
        !               we assume that pointwise data follows a Gaussian Process with a covariance function
        !               defined by the squared exponential (SE) kernel k(a,b). In other words, this
        !               function returns the integral of k(a,p2) with respect to a
        !               over the volume [x1-dx, x1+dx]x[y1-dy,y1+dy].
        !               
        !               See Reyes et al. 2018 "A New Class of High-Order Methods for Fluid Dynamics
        !               Simulations Using Gaussian Process Modeling: One-Dimensional Case" equation
        !               equation 23 or Bourgeois and Lee 2022 "GP-MOOD: A positivity-preserving high-order
        !               finite volume method for hyperbolic conservation laws" equations 10 and 11.
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
            r1 = (delta+0.5_qp)/(SQRT(2)*ell_over_dl)
            r2 = (delta-0.5_qp)/(SQRT(2)*ell_over_dl)
            kernel = kernel*&
                SQRT(qp_pi/2._qp) * (ell/dl(i)) * &
                (ERF(r1) - ERF(r2))
        end do

    end function GP_volAvgToPointSqExpKernel_

end module GP_volAvgToPointSqExpKernel
