module GP_volAvgToPointPredVect
    
    use definitions, only: qp, ndim
    use GP_volAvgToVolAvgSqExpKernel
    use GP_volAvgToPointSqExpKernel

    implicit none

contains

    function GP_volAvgToPointPredVect_(XX, xstr, nStenc, nPred, dl, ell) result(predVect)
        ! function:     GP_volAvgToPointPredVect_
        ! purpose:      
        !               
        !               
        ! 
        ! Inputs:       
        !               
        ! Outputs:      
        ! ------------------------------------------------------------

        implicit none
        integer, intent(in) :: nStenc, nPred
        real(qp), intent(in), dimension(ndim, nStenc) :: XX
        real(qp), intent(in), dimension(ndim, nPred) :: xstr
        real(qp), intent(in) :: dl(ndim), ell(ndim)
        real(qp) :: predVect(nPred, nStenc)

        ! local variables
        integer :: i, j
        real(qp) :: Kcov(nStenc, nStenc), Kstr(nPred, nStenc)

        do i=1,nStenc
            do j=1,nStenc
                Kcov(i,j) = GP_volAvgToVolAvgSqExpKernel_(XX(:,i), XX(:,j), dl, ell)
            end do
        end do

        do i=1,nPred
            do j=1,nStenc
                Kstr(i,j) = GP_volAvgToPointSqExpKernel_(XX(:,i), XX(:,j), dl, ell)
            end do
        end do

    end function GP_volAvgToPointPredVect_

end module GP_volAvgToPointPredVect
