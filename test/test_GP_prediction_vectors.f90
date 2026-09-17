program test_GP_prediction_vectors
    ! program:      test_GP_prediction_vectors

    use definitions, only: max_string_length
    use GP, only: GP_predictionVectors, GP_init, GP_finalize, GP_maxRadius, &
        GP_nStenc, GP_nStencMax, GP_nPred, GP_nPredMax
    use assert, only: assert_close, assert_summary

    implicit none

    character(len=max_string_length) :: paramfile
    real, allocatable :: oneNorm(:)
    real :: expected_sum_of_oneNorm
    integer :: exitStat
    integer :: rr

    paramfile = "par/test_GP_prediction_vectors.par"
    call GP_init(paramfile)

    allocate(oneNorm(GP_nPredMax))
    
    do rr=1,GP_maxRadius
        expected_sum_of_oneNorm = GP_nPred(rr)
        oneNorm = SUM(GP_predictionVectors(:,:,rr), DIM=1)
        call assert_close(SUM(oneNorm), expected_sum_of_oneNorm, 1e-5, "GP_prediction vectors all have 1-norm equal to unity")
    end do

    call assert_summary(exitStat)
    call EXIT(exitStat)

end program test_GP_prediction_vectors
