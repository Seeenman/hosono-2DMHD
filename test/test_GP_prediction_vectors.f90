program test_GP_prediction_vectors
    ! program:      test_GP_prediction_vectors

    use definitions, only: max_string_length
    use GP, only: GP_predictionVectors, GP_init, GP_finalize, GP_maxRadius, &
        GP_nStenc, GP_nStencMax, GP_nPred, GP_nPredMax
    use grid, only: grid_init, grid_finalize
    use assert, only: assert_close, assert_summary

    implicit none

    character(len=max_string_length) :: paramfile
    real :: oneNorm
    integer :: exitStat
    integer :: rr, i
    real, allocatable :: vec(:)

    paramfile = "par/test_GP_prediction_vectors.par"
    call GP_init(paramfile)

    do rr=1,GP_maxRadius
        do i=1,GP_nPred(rr)
            oneNorm = SUM(GP_predictionVectors(:,i,rr))
            call assert_close(oneNorm, 1.0, 1e-2, "GP_prediction vector has 1-norm close to unity")
        end do
    end do

    ! check that 
    allocate(vec(GP_nStencMax))

    deallocate(vec)

    call GP_finalize()

    call assert_summary(exitStat)
    call EXIT(exitStat)

end program test_GP_prediction_vectors
