module GP_volAvgToPointPredVect
    
    use definitions, only: qp, ndim
    use GP_volAvgToVolAvgSqExpKernel
    use GP_volAvgToPointSqExpKernel
    use linAlgQuadPrecision, only: linAlgQuadPrecision_solveSPD

    implicit none

contains

    function GP_volAvgToPointPredVect_(XX, XXstr, nStenc, nPred, dl, ell) result(predVect)
        ! function:     GP_volAvgToPointPredVect_
        ! purpose:      Construct the Gaussian process prediction vector(s) for predicting 
        !               pointwise data given volume averaged data.
        !               
        ! Inputs:       - XX(ndim, nStenc) (quadruple precision real) the cell
        !                 center coordinates of all of the cells in the GP stencil
        !               - XXstr(ndim, nPred) (quadruple precision real) the coordinates of the points
        !                 at which we want to predict pointwise data
        !               - nStenc (integer) the number of cells in the GP stencil
        !               - nPred (integer) the number of points at which we want to predict pointwise data
        !               - dl(ndim) (quadruple precision real) the grid spacing. dl(1) = dx, dl(2) = dy
        !               - ell is the hyperparameter of the squared exponential kernel
        !               
        ! Outputs:      - predVect(nStenc, nPred) the GP prediction vector(s).
        !                 The ith column of predVect is the prediction vector for
        !                 the ith point at which we want to predict pointwise data.
        ! ------------------------------------------------------------

        implicit none
        integer, intent(in) :: nStenc, nPred
        real(qp), intent(in), dimension(ndim, nStenc) :: XX
        real(qp), intent(in), dimension(ndim, nPred) :: XXstr
        real(qp), intent(in) :: dl(ndim), ell(ndim)
        real(qp) :: predVect(nStenc, nPred)

        ! local variables
        integer :: i, j
        real(qp) :: Kcov(nStenc, nStenc), Kstr(nStenc, nPred)

        do i=1,nStenc
            do j=1,nStenc
                Kcov(i,j) = GP_volAvgToVolAvgSqExpKernel_(XX(:,i), XX(:,j), dl, ell)
            end do
        end do

        do i=1,nStenc
            do j=1,nPred
                Kstr(i,j) = GP_volAvgToPointSqExpKernel_(XX(:,i), XXstr(:,j), dl, ell)
            end do
        end do

        predVect = linAlgQuadPrecision_solveSPD(Kcov, Kstr)

    end function GP_volAvgToPointPredVect_

end module GP_volAvgToPointPredVect
