module GP

    use definitions, only: ndim, qp, qp_pi, max_string_length, xdir, ydir
    use readParamFile, only: readParamFile_int, readParamFile_quadPrecisionReal
    use linAlgQuadPrecision, only: linAlgQuadPrecision_solveSPD

    implicit none

    private

    ! public variables
    public :: GP_maxRadius
    public :: GP_nPred
    public :: GP_nStenc
    public :: GP_nQuadrature
    public :: GP_nPredMax
    public :: GP_nStencMax
    public :: GP_nQuadratureMax
    public :: GP_predictionVectors
    public :: GP_quadratureWeights
    public :: GP_stencIdxs

    ! public routines
    public :: GP_init
    public :: GP_finalize

    ! Gaussian Process Global variables
    integer :: GP_maxRadius
    integer:: GP_nStencMax, GP_nPredMax, GP_nQuadratureMax
    integer, allocatable, dimension(:) :: GP_nStenc, GP_nPred, GP_nQuadrature
    integer, allocatable :: GP_stencIdxs(:,:,:)
    real, allocatable :: GP_predictionVectors(:,:,:)
    real, allocatable, dimension(:,:) :: GP_quadraturePoints, GP_quadratureWeights

contains

    subroutine GP_init(paramfile)
        ! subroutine:   GP_init
        ! purpose:      
        ! 
        ! Inputs:       
        !               
        ! Outputs:      
        ! ------------------------------------------------------------
        implicit none
        character(len=max_string_length), intent(in) :: paramfile
        ! local variables
        real(qp), allocatable :: XX(:,:), XXstr(:,:), predVect(:,:), quadraturePoints(:)
        real(qp), dimension(ndim) :: domainBeg, domainEnd, nGrid, dl
        real(qp) :: ell
        integer :: i, j, rr, counter

        write(*,*) "=============================================================="
        write(*,*) "Initializing Gaussian Process variables:"
        write(*,*) "--------------------------------------------------------------"

        domainBeg(xdir) = readParamFile_quadPrecisionReal(paramfile, "grid_xBeg")
        domainEnd(xdir) = readParamFile_quadPrecisionReal(paramfile, "grid_xEnd")
        domainBeg(ydir) = readParamFile_quadPrecisionReal(paramfile, "grid_yBeg")
        domainEnd(ydir) = readParamFile_quadPrecisionReal(paramfile, "grid_yEnd")
        nGrid(xdir) = readParamFile_int(paramfile, "grid_Nx")
        nGrid(ydir) = readParamFile_int(paramfile, "grid_Ny")
        dl(xdir) = (domainEnd(xdir)-domainBeg(xdir))/nGrid(xdir)
        dl(ydir) = (domainEnd(ydir)-domainBeg(ydir))/nGrid(ydir)
        ell = readParamFile_quadPrecisionReal(paramfile, "GP_ellOverDelta")*MINVAL(dl)
        print*, "ell: ", ell

        GP_maxRadius = readParamFile_int(paramfile, "GP_maxRadius")
        allocate(GP_nStenc(GP_maxRadius))
        allocate(GP_nQuadrature(GP_maxRadius))
        allocate(GP_nPred(GP_maxRadius))
        do rr = 1,GP_maxRadius
            GP_nStenc(rr) = 1+2*rr*(rr+1) ! number of cells in the GP stencil at each order
            GP_nQuadrature(rr) = rr+1 ! number of quadrature points per face
            GP_nPred(rr) = 4*GP_nQuadrature(rr) ! number of points at which we need to predict (all of the quadrature points)
        end do
        GP_nStencMax = GP_nStenc(GP_maxRadius)
        GP_nQuadratureMax = GP_nQuadrature(GP_maxRadius)
        GP_nPredMax = GP_nPred(GP_maxRadius)

        ! fill in GP stencil indices
        allocate(GP_stencIdxs(ndim, GP_nStencMax, GP_maxRadius))
        GP_stencIdxs = 0
        do rr=1,GP_maxRadius
            counter=0
            do i=-rr, rr
                do j=-rr, rr
                    if (ABS(i)+ABS(j)<=rr) then
                        counter = counter + 1
                        GP_stencIdxs(:,counter,rr) = [i,j]
                    end if
                end do
            end do
        end do

        ! Allocate quadrature coordinates and quadrature weights
        ! based on GP_maxRadius
        allocate(GP_quadraturePoints(GP_nQuadratureMax, GP_maxRadius))
        allocate(GP_quadratureWeights(GP_nQuadratureMax, GP_maxRadius))
        allocate(GP_predictionVectors(GP_nStencMax, GP_nPredMax, GP_maxRadius))
        GP_quadraturePoints = 0.0 !TODO this might let some bugs through silently

        do rr=1,GP_maxRadius
            ! fill in GP training inputs
            allocate(XX(ndim, GP_nStenc(rr)))
            allocate(XXstr(ndim, GP_nPred(rr)))
            allocate(quadraturePoints(GP_nQuadrature(rr)))
            XX(xdir, :) = GP_stencIdxs(xdir, 1:GP_nStenc(rr),rr)*dl(xdir)
            XX(ydir, :) = GP_stencIdxs(ydir, 1:GP_nStenc(rr),rr)*dl(ydir)

            ! set quadrature coordinates and quadrature weights
            ! NOTE: quadarature weights don't have to be quadruple precision
            if (rr == 1) then
                ! GP radius of 1
                ! GP spatial order of accuracy = 2*1+1 = 3
                ! use 4th order, 2 point quadrature rule
                quadraturePoints(1) = 1.0_qp/2.0_qp/SQRT(3.0_qp)
                quadraturePoints(2) = -quadraturePoints(1)
                GP_quadratureWeights(1, rr) = 1.0/2.0
                GP_quadratureWeights(2, rr) = GP_quadratureWeights(1, rr)
            else if (rr == 2) then
                ! GP radius of 2
                ! GP spatial order of accuracy = 2*2+1 = 5
                ! use 6th order, 3 point quadrature rule
                quadraturePoints(1) = 1.0_qp/2.0_qp*SQRT(3.0_qp/5.0_qp)
                quadraturePoints(2) = 0.0_qp
                quadraturePoints(3) = -quadraturePoints(1)
                GP_quadratureWeights(1, rr) = 5.0/18.0
                GP_quadratureWeights(2, rr) = 8.0/18.0
                GP_quadratureWeights(3, rr) = GP_quadratureWeights(1, rr)
            else if (rr == 3) then
                ! GP radius of 3
                ! GP spatial order of accuracy = 2*3+1 = 7
                ! use 8th order, 4 point quadrature rule
                quadraturePoints(1) = 1.0_qp/2.0_qp*SQRT(3.0_qp/7.0_qp+2.0_qp/7.0_qp*SQRT(6.0_qp/5.0_qp))
                quadraturePoints(2) = 1.0_qp/2.0_qp*SQRT(3.0_qp/7.0_qp-2.0_qp/7.0_qp*SQRT(6.0_qp/5.0_qp))
                quadraturePoints(3) = -quadraturePoints(2)
                quadraturePoints(4) = -quadraturePoints(1)
                GP_quadratureWeights(1, rr) = (18.0-SQRT(30.0))/72.0
                GP_quadratureWeights(2, rr) = (18.0+SQRT(30.0))/72.0
                GP_quadratureWeights(3, rr) = GP_quadratureWeights(2, rr)
                GP_quadratureWeights(4, rr) = GP_quadratureWeights(1, rr)
            else
                error stop "GP radius larger than 3 not currently supported"
            end if

            ! fill in GP test outputs based on quadrature points
            ! upper face
            XXstr(xdir, 0*GP_nQuadrature(rr)+1:1*GP_nQuadrature(rr)) = quadraturePoints*dl(xdir)
            XXstr(ydir, 0*GP_nQuadrature(rr)+1:1*GP_nQuadrature(rr)) = dl(ydir)/2.0_qp
            ! lower face
            XXstr(xdir, 1*GP_nQuadrature(rr)+1:2*GP_nQuadrature(rr)) = quadraturePoints*dl(xdir)
            XXstr(ydir, 1*GP_nQuadrature(rr)+1:2*GP_nQuadrature(rr)) = -dl(ydir)/2.0_qp
            ! right face
            XXstr(xdir, 2*GP_nQuadrature(rr)+1:3*GP_nQuadrature(rr)) = dl(xdir)/2.0_qp
            XXstr(ydir, 2*GP_nQuadrature(rr)+1:3*GP_nQuadrature(rr)) = quadraturePoints*dl(ydir)
            ! left face
            XXstr(xdir, 3*GP_nQuadrature(rr)+1:4*GP_nQuadrature(rr)) = -dl(xdir)/2.0_qp
            XXstr(ydir, 3*GP_nQuadrature(rr)+1:4*GP_nQuadrature(rr)) = quadraturePoints*dl(ydir)

            allocate(predVect(GP_nStenc(rr), GP_nPred(rr)))
            predVect = GP_volAvgToPointPredVect(XX, XXstr, GP_nStenc(rr), GP_nPred(rr), dl, ell)

            make 1-norm of each individual prediction vector equal to 1
            do i=1,GP_nPred(rr)
                predVect(:,i) = predVect(:,i)/SUM(predVect(:,i))
            end do

            ! convert from quad precision to double precision
            GP_predictionVectors(1:GP_nStenc(rr), 1:GP_nPred(rr), rr) = REAL(predVect) 
            GP_quadraturePoints(1:GP_nQuadrature(rr), rr) = REAL(quadraturePoints)

            deallocate(quadraturePoints)
            deallocate(XX)
            deallocate(XXstr)
            deallocate(predVect)
        end do

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "GP variables initialized"
        write(*,*) "=============================================================="

    end subroutine GP_init

    subroutine GP_finalize()
        implicit none
        deallocate(GP_predictionVectors)
        deallocate(GP_stencIdxs)
        deallocate(GP_quadraturePoints)
        deallocate(GP_quadratureWeights)
        deallocate(GP_nStenc)
        deallocate(GP_nQuadrature)
        deallocate(GP_nPred)
        write(*,*) "=============================================================="
        write(*,*) "GP variables deallocated."
        write(*,*) "=============================================================="
    end subroutine GP_finalize

    pure function GP_volAvgToVolAvgSqExpKernel(p1, p2, dl, ell) result(kernel)
        ! function:     GP_volAvgToVolAvgSqExpKernel
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
        real(qp) :: sep
        real(qp) :: kernel
        real(qp) :: sqrt2_ell, two_ell_sq, sqrt_pi_ell_sq
        real(qp) :: r1, r2, r3, r4, r5, r6
        integer :: i
        real(qp), parameter :: sqrt_two = SQRT(2._qp)

        sqrt2_ell = sqrt_two*ell
        two_ell_sq = 2._qp*ell**2
        sqrt_pi_ell_sq = SQRT(qp_pi)*ell**2

        kernel = 1._qp
        do i=1,ndim
            sep = p1(i)-p2(i) ! numerator of equation 13 in Bourgeois and Lee
            r1 = (sep+dl(i))/sqrt2_ell
            r2 = (sep-dl(i))/sqrt2_ell
            r3 = - (sep+dl(i))**2 / two_ell_sq
            r4 = - (sep-dl(i))**2 / two_ell_sq
            r5 = sep/sqrt2_ell
            r6 = - sep**2 / two_ell_sq

            kernel = kernel*&
                sqrt_pi_ell_sq/dl(i)**2 * &
                (&
                r1*ERF(r1) + r2*ERF(r2) &
                + (EXP(r3) + EXP(r4))/SQRT(qp_pi) &
                - 2._qp*(r5*ERF(r5) + EXP(r6)/SQRT(qp_pi)) &
                )
        end do

    end function GP_volAvgToVolAvgSqExpKernel

    pure function GP_volAvgToPointSqExpKernel(p1, p2, dl, ell) result(kernel)
        ! function:     GP_volAvgToPointSqExpKernel
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
        real(qp) :: sep
        real(qp) :: kernel
        real(qp) :: sqrt2_ell, half_delta
        real(qp) :: r1, r2
        integer :: i
        real(qp), parameter :: sqrt_two     = SQRT(2._qp)
        real(qp), parameter :: sqrt_half_pi = SQRT(qp_pi/2._qp)


        sqrt2_ell = sqrt_two*ell

        kernel = 1._qp
        do i=1,ndim
            half_delta = dl(i)/2._qp
            sep = p1(i)-p2(i) ! numerator of equation 13 in Bourgeois and Lee
            r1 = (sep+half_delta)/sqrt2_ell
            r2 = (sep-half_delta)/sqrt2_ell
            kernel = kernel*&
                sqrt_half_pi * (ell/dl(i)) * (ERF(r1) - ERF(r2))
        end do

    end function GP_volAvgToPointSqExpKernel

    function GP_volAvgToPointPredVect(XX, XXstr, nStenc, nPred, dl, ell) result(predVect)
        ! function:     GP_volAvgToPointPredVect
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
        real(qp), intent(in) :: dl(ndim), ell
        real(qp) :: predVect(nStenc, nPred)

        ! local variables
        integer :: i, j
        real(qp) :: Kcov(nStenc, nStenc), Kstr(nStenc, nPred)

        do i=1,nStenc
            do j=1,nStenc
                Kcov(i,j) = GP_volAvgToVolAvgSqExpKernel(XX(:,i), XX(:,j), dl, ell)
            end do
        end do

        do i=1,nStenc
            do j=1,nPred
                Kstr(i,j) = GP_volAvgToPointSqExpKernel(XX(:,i), XXstr(:,j), dl, ell)
            end do
        end do

        predVect = linAlgQuadPrecision_solveSPD(Kcov, Kstr)

    end function GP_volAvgToPointPredVect

end module GP
