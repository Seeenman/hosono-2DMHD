module linAlgQuadPrecision

    use definitions, only: qp

    implicit none

    interface linAlgQuadPrecision_forwardSub
        module procedure linAlgQuadPrecision_forwardSub_matRHS
        module procedure linAlgQuadPrecision_forwardSub_vecRHS
    end interface linAlgQuadPrecision_forwardSub

    interface linAlgQuadPrecision_backSub
        module procedure linAlgQuadPrecision_backSub_matRHS
        module procedure linAlgQuadPrecision_backSub_vecRHS
    end interface linAlgQuadPrecision_backSub

    interface linAlgQuadPrecision_solve
        module procedure linAlgQuadPrecision_solve_matRHS
        module procedure linAlgQuadPrecision_solve_vecRHS
    end interface linAlgQuadPrecision_solve

    interface linAlgQuadPrecision_solveSPD
        module procedure linAlgQuadPrecision_solveSPD_matRHS
        module procedure linAlgQuadPrecision_solveSPD_vecRHS
    end interface linAlgQuadPrecision_solveSPD

contains

    subroutine linAlgQuadPrecision_LU(A, L, U, p, m)
        ! LU factorization of A via Gaussian elimination with partial pivoting.
        ! Given an mxm input matrix A, returns the factorization
        ! PA = LU where
        ! L is a lower-triangular matrix
        ! U is an upper-triangular matrix
        ! P is a permutation matrix
        ! 
        ! note that instead of the entire matrix P being returned directly
        ! we instead return a permutation vector p such that
        ! p(i) = position of the 1 in the ith row of P
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in) :: A(m,m)
        real(qp), intent(out) :: L(m,m), U(m,m)
        integer, intent(out) :: p(m)

        integer :: i,j,k,s
        real(qp) :: v(SIZE(A,1))

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

    function linAlgQuadPrecision_choleskyDecomp(A,m) result(R)
        ! cholesky decomposition of an mxm real symmetric
        ! positive definite matrix
        !
        ! returns R in the upper-triangle of A
        ! where the cholesky decomposition of A is
        ! A = R'R 
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in) :: A(m,m)
        real(qp) :: R(m,m)
        integer :: k,j
        
        R = A
        do k=1,m
            do j=k+1,m
                R(j,j:m) = R(j,j:m) - R(k,j:m)*R(k,j)/R(k,k)
            end do
            R(k,k:m) = R(k,k:m)/SQRT(R(k,k))
        end do

    end function linAlgQuadPrecision_choleskyDecomp

    function linAlgQuadPrecision_forwardSub_matRHS(L, m, b, n) result(y)
        ! Solves an mxm lower-triangular system Ly=b for y.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ly=b.
        ! The result y has the same shape as b.
        implicit none
        integer, intent(in) :: m, n
        real(qp), intent(in) :: L(m,m), b(m,n)
        real(qp) :: y(m,n)
        real(qp) :: ss(n)
        integer :: i,j
        
        y(1,:) = b(1,:)/L(1,1)
        do i=2,m
            ss=0._qp
            do j=1,i-1
                ss = ss + L(i,j)*y(j,:)
            end do
            y(i,:) = (b(i,:) - ss)/L(i,i)
        end do
    end function linAlgQuadPrecision_forwardSub_matRHS

    function linAlgQuadPrecision_forwardSub_vecRHS(L, m, b) result(y)
        ! Solves an mxm lower-triangular system Ly=b for y.
        ! b is an mx1 column vector
        ! The result y has the same shape as b.
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in) :: L(m,m), b(m)
        real(qp) :: y(m)
        real(qp) :: ss
        integer :: i,j
        
        y(1) = b(1)/L(1,1)
        do i=2,m
            ss=0._qp
            do j=1,i-1
                ss = ss + L(i,j)*y(j)
            end do
            y(i) = (b(i) - ss)/L(i,i)
        end do
    end function linAlgQuadPrecision_forwardSub_vecRHS

    function linAlgQuadPrecision_backSub_matRHS(U, m, y, n) result(x)
        ! Solves an mxm upper-triangular system Ux=y for x
        ! y is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ux=y.
        ! The result x has the same shape as y.
        implicit none
        integer, intent(in) :: m, n
        real(qp), intent(in) :: U(m,m), y(m,n)
        real(qp) :: x(m,n)
        real(qp) :: ss(n)
        integer :: i,j
        
        x(m,:) = y(m,:)/U(m,m)
        do i=m-1,1,-1
            ss=0._qp
            do j=i+1,m
                ss = ss + U(i,j)*x(j,:)
            end do
            x(i,:) = (y(i,:) - ss)/U(i,i)
        end do
    end function linAlgQuadPrecision_backSub_matRHS

    function linAlgQuadPrecision_backSub_vecRHS(U, m, y) result(x)
        ! Solves an mxm upper-triangular system Ux=y for x
        ! y is an mx1 column vector
        ! The result x has the same shape as y.
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in) :: U(m,m), y(m)
        real(qp) :: x(m)
        real(qp) :: ss
        integer :: i,j
        
        x(m) = y(m)/U(m,m)
        do i=m-1,1,-1
            ss=0._qp
            do j=i+1,m
                ss = ss + U(i,j)*x(j)
            end do
            x(i) = (y(i) - ss)/U(i,i)
        end do
    end function linAlgQuadPrecision_backSub_vecRHS

    function linAlgQuadPrecision_solve_matRHS(A, m, b, n) result(x)
        ! Solves the mxm system Ax=b for x.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ax=b.
        ! The result x has the same shape as b.
        implicit none
        integer, intent(in) :: m, n
        real(qp), intent(in)  :: A(m,m), b(m,n)
        real(qp) :: x(m,n)
        real(qp) :: L(m,m), U(m,m)
        integer :: p(m)
        real(qp) :: Pb(m,n)
        integer :: i

        call linAlgQuadPrecision_LU(A, L, U, p, m)

        ! apply permuation vector to b
        do i=1,m
            Pb(i,:) = b(p(i), :)
        end do

        ! now we are set up to solve LUx = PAx = Pb for x.
        x = linAlgQuadPrecision_forwardSub(L, m, Pb, n) ! first solve Lx' = Pb where x'=Ux
        x = linAlgQuadPrecision_backSub(U, m, x, n) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solve_matRHS

    function linAlgQuadPrecision_solve_vecRHS(A, m, b) result(x)
        ! Solves the mxm system Ax=b for x.
        ! b is an mx1 column vector
        ! The result x has the same shape as b.
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in)  :: A(m,m), b(m)
        real(qp) :: x(m)
        real(qp) :: L(m,m), U(m,m)
        integer :: p(m)
        real(qp) :: Pb(m)
        integer :: i

        call linAlgQuadPrecision_LU(A, L, U, p, m)

        ! apply permuation vector to b
        do i=1,m
            Pb(i) = b(p(i))
        end do

        ! now we are set up to solve LUx = PAx = Pb for x.
        x = linAlgQuadPrecision_forwardSub(L, m, Pb) ! first solve Lx' = Pb where x'=Ux
        x = linAlgQuadPrecision_backSub(U, m, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solve_vecRHS

    function linAlgQuadPrecision_solveSPD_matRHS(A, m, b, n) result(x)
        ! Solves the system Ax=b for x when A is symmetric positive definite.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ax=b.
        ! The result x has the same shape as b.
        implicit none
        integer, intent(in) :: m, n
        real(qp), intent(in) :: A(m,m)
        real(qp), intent(in) :: b(m,n)
        real(qp) :: x(m,n)
        real(qp) :: R(m,m), Rt(m,m)

        R = linAlgQuadPrecision_choleskyDecomp(A, m)
        Rt = TRANSPOSE(R)

        ! now we are set up to solve LUx = Ax = b for x where L=R^T and U=R
        x = linAlgQuadPrecision_forwardSub(Rt, m, b, n) ! first solve Lx' = b where x'=Ux
        x = linAlgQuadPrecision_backSub(R, m, x, n) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solveSPD_matRHS

    function linAlgQuadPrecision_solveSPD_vecRHS(A, m, b) result(x)
        ! Solves the system Ax=b for x when A is symmetric positive definite.
        ! b is an mx1 column vector
        ! The result x has the same shape as b.
        implicit none
        integer, intent(in) :: m
        real(qp), intent(in) :: A(m,m)
        real(qp), intent(in) :: b(m)
        real(qp) :: x(m)
        real(qp) :: R(m,m), Rt(m,m)

        R = linAlgQuadPrecision_choleskyDecomp(A, m)
        Rt = TRANSPOSE(R)

        ! now we are set up to solve LUx = Ax = b for x where L=R^T and U=R
        x = linAlgQuadPrecision_forwardSub(Rt, m, b) ! first solve Lx' = b where x'=Ux
        x = linAlgQuadPrecision_backSub(R, m, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solveSPD_vecRHS

end module linAlgQuadPrecision
