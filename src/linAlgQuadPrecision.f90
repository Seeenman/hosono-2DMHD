module linAlgQuadPrecision

    use definitions, only: qp

    implicit none

contains

    function linAlgQuadPrecision_choleskyDecomp(A) result(R)
        ! cholesky decomposition of a real symmetric
        ! positive definite matrix
        !
        ! returns R in the upper-triangle of A
        ! where the cholesky decomposition of A is
        ! A = R'R 
        implicit none
        real(qp), intent(in)  :: A(:,:)
        real(qp) :: R(SIZE(A,1),SIZE(A,1))
        integer :: k,j,m
        
        m = SIZE(A,1)
        R = A
        do k=1,m
            do j=k+1,m
                R(j,j:m) = R(j,j:m) - R(k,j:m)*R(k,j)/R(k,k)
            end do
            R(k,k:m) = R(k,k:m)/SQRT(R(k,k))
        end do

    end function linAlgQuadPrecision_choleskyDecomp

    function linAlgQuadPrecision_forwardSub(L, b) result(y)
        ! Solves a lower-triangular system Ly=b for y
        implicit none
        real(qp), intent(in)  :: L(:,:), b(:)
        real(qp) :: y(SIZE(L,1))
        real(qp) :: ss
        integer :: i,j,m
        
        m = SIZE(L,1)
        y(1) = b(1)/L(1,1)
        do i=2,m
            ss=0._qp
            do j=1,i-1
                ss = ss + L(i,j)*y(j)
            end do
            y(i) = (b(i) - ss)/L(i,i)
        end do
    end function linAlgQuadPrecision_forwardSub

    function linAlgQuadPrecision_backSub(U, y) result(x)
        ! Solves an upper-triangular system Ux=y for x
        implicit none
        real(qp), intent(in)  :: U(:,:), y(:)
        real(qp) :: x(SIZE(U,1))
        real(qp) :: ss
        integer :: i,j,m
        
        m = SIZE(U,1)
        x(m) = y(m)/U(m,m)
        do i=m-1,1,-1
            ss=0._qp
            do j=i+1,m
                ss = ss + U(i,j)*x(j)
            end do
            x(i) = (y(i) - ss)/U(i,i)
        end do
    end function linAlgQuadPrecision_backSub

    subroutine linAlgQuadPrecision_LU(A, L, U, p)
        ! LU factorization of A via Gaussian elimination with partial pivoting.
        ! Given input matrix A, returns the factorization
        ! PA = LU where
        ! L is a lower-triangular matrix
        ! U is an upper-triangular matrix
        ! P is a permutation matrix
        ! 
        ! note that instead of the entire matrix P being returned directly
        ! we instead return a permutation vector p such that
        ! p(i) = position of the 1 in the ith row of P
        implicit none
        real(qp), intent(in)  :: A(:,:)
        real(qp), intent(out)  :: L(SIZE(A,1),SIZE(A,1)), U(SIZE(A,1),SIZE(A,1))
        integer, intent(out) :: p(SIZE(A,1))

        integer :: i,j,k,m,s
        real(qp) :: v(SIZE(A,1))

        m = SIZE(A,1)

        ! check that A is square
        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_LU: A must be square"
        end if

        U = A
        L = 0._qp
        p = 0
        do i=1,m
            L(i,i) = 1._qp
            p(i) = i
        end do

        do k=1,m-1
            i = MAXLOC(ABS(U(k:,k)), DIM=1) + k-1

            if (i>k) then
                ! interchange rows
                v(k:) = U(k,k:)
                U(k,k:) = U(i,k:)
                U(i,k:) = v(k:)

                v(1:k-1) = L(k,1:k-1)
                L(k,1:k-1) = L(i,1:k-1)
                L(i,1:k-1) = v(1:k-1)

                s = p(k)
                p(k) = p(i)
                p(i) = s
            end if
            do j=k+1, m
                L(j,k) = U(j,k)/U(k,k)
                U(j,k+1:) = U(j,k+1:) - L(j,k)*U(k,k+1:)
                U(j,k) = 0._qp
            end do
        end do
        
    end subroutine linAlgQuadPrecision_LU

    function linAlgQuadPrecision_solveLinearSystem(A, b) result(x)
        ! Solves the system Ax=b for x.
        ! If A is not invertible then prints an error message and stops
        implicit none
        real(qp), intent(in)  :: A(:,:), b(:)
        real(qp) :: x(SIZE(A,1))
        real(qp) :: L(SIZE(A,1),SIZE(A,1)), U(SIZE(A,1),SIZE(A,1))
        integer :: p(SIZE(A,1))
        real(qp) :: det, Pb(SIZE(A,1))
        integer :: m, i

        m = size(A,1)

        call linAlgQuadPrecision_LU(A, L, U, p)

        ! apply permuation vector to b
        do i=1,m
            Pb(i) = b(p(i))
        end do

        ! now we are set up to solve LUx = PAx = Pb for x.
        x = linAlgQuadPrecision_forwardSub(L, Pb) ! first solve Lx' = Pb where x'=Ux
        x = linAlgQuadPrecision_backSub(U, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solveLinearSystem

end module linAlgQuadPrecision
