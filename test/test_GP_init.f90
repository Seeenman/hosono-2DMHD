program test_GP_init
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use definitions, only: ndim, max_string_length
    use GP, only: GP_init, GP_finalize, GP_maxRadius, GP_nStenc, GP_nQuadrature, GP_nPred, &
                  GP_nStencMax, GP_nQuadratureMax, GP_nPredMax, &
                  GP_stencIdxs, GP_predictionVectors, GP_quadratureWeights
    use assert, only: assert_true, assert_equal_int, assert_close, assert_summary
    implicit none

    character(len=max_string_length) :: paramfile
    integer :: exitStat
    integer :: expectedRadius, expectedNStenc, expectedNQuad, expectedNPred
    real :: weightSum

    write(*,*) "================================================================"
    write(*,*) "test_GP_init"
    write(*,*) "================================================================"

    paramfile = "par/test_GP_init.par"
    call GP_init(paramfile)

    ! values matching test/par/test_GP_init.par
    expectedRadius = 2
    expectedNQuad  = expectedRadius + 1
    expectedNStenc = 1 + 2*expectedRadius*(expectedRadius+1)
    expectedNPred  = 4*expectedNQuad

    call assert_equal_int(GP_maxRadius, expectedRadius, "GP_maxRadius matches param file")
    call assert_equal_int(GP_nQuadratureMax, expectedNQuad, "GP_nQuadratureMax = GP_maxRadius+1")
    call assert_equal_int(GP_nStencMax, expectedNStenc, "GP_nStencMax = 1+2*R*(R+1)")
    call assert_equal_int(GP_nPredMax, expectedNPred, "GP_nPredMax = 4*GP_nQuadrature")

    call assert_true(ALLOCATED(GP_stencIdxs), "GP_stencIdxs is allocated")
    call assert_equal_int(SIZE(GP_stencIdxs,1), ndim, "GP_stencIdxs first dim = ndim")
    call assert_equal_int(SIZE(GP_stencIdxs,2), expectedNStenc, "GP_stencIdxs second dim = GP_nStenc")
    call assert_equal_int(COUNT(GP_stencIdxs(1,:, expectedRadius)==0 .and. GP_stencIdxs(2,:, expectedRadius)==0), 1, &
        "GP_stencIdxs contains the center cell (0,0) exactly once")

    call assert_true(ALLOCATED(GP_predictionVectors), "GP_predictionVectors is allocated")
    call assert_equal_int(SIZE(GP_predictionVectors,1), expectedNStenc, "GP_predictionVectors rows = GP_nStenc")
    call assert_equal_int(SIZE(GP_predictionVectors,2), expectedNPred, "GP_predictionVectors cols = GP_nPred")
    call assert_true(.not. ANY(ieee_is_nan(GP_predictionVectors)), &
        "GP_predictionVectors contains no NaNs")

    call assert_true(ALLOCATED(GP_quadratureWeights), "GP_quadratureWeights is allocated")
    call assert_equal_int(SIZE(GP_quadratureWeights(:,GP_maxRadius)), expectedNQuad, "GP_quadratureWeights size = GP_nQuadrature")
    weightSum = SUM(GP_quadratureWeights(:,GP_maxRadius))
    call assert_close(weightSum, 1.0, 1.0e-5, "GP_quadratureWeights sum to 1")

    call GP_finalize()

    call assert_summary(exitStat)
    call EXIT(exitStat)

end program test_GP_init
