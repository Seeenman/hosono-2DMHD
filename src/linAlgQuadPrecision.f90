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

    subroutine linAlgQuadPrecision_LU(A, L, U, p)
        ! LU factorization of A via Gaussian elimination with partial pivoting.
        ! Given a square input matrix A, returns the factorization
        ! PA = LU where
        ! L is a lower-triangular matrix
        ! U is an upper-triangular matrix
        ! P is a permutation matrix
        !
        ! note that instead of the entire matrix P being returned directly
        ! we instead return a permutation vector p such that
        ! p(i) = position of the 1 in the ith row of P
        !
        ! L, U and p are intent(out) assumed-shape arguments, so the caller
        ! must supply them already sized (m,m), (m,m) and (m) respectively.
        implicit none
        real(qp), intent(in)  :: A(:,:)
        real(qp), intent(out) :: L(:,:), U(:,:)
        integer,  intent(out) :: p(:)

        integer  :: i, j, k, s, m
        real(qp) :: v(SIZE(A,1))

        m = SIZE(A,1)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_LU: A must be square"
        end if
        if (ANY(SHAPE(L) /= [m,m])) then
            error stop "linAlgQuadPrecision_LU: L must have the same shape as A"
        end if
        if (ANY(SHAPE(U) /= [m,m])) then
            error stop "linAlgQuadPrecision_LU: U must have the same shape as A"
        end if
        if (SIZE(p) /= m) then
            error stop "linAlgQuadPrecision_LU: p must have length SIZE(A,1)"
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

    function linAlgQuadPrecision_choleskyDecomp(A) result(R)
        ! Cholesky decomposition of a real symmetric positive definite matrix A.
        !
        ! Returns upper-triangular R such that A = R^T R.
        !
        ! Only the upper triangle of A is read; the strict lower triangle of the
        ! input is ignored, and the strict lower triangle of R is set to zero.
        implicit none
        real(qp), intent(in) :: A(:,:)
        real(qp) :: R(SIZE(A,1), SIZE(A,1))

        integer :: k, j, m

        m = SIZE(A,1)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_choleskyDecomp: A must be square"
        end if
        if (m < 1) then
            error stop "linAlgQuadPrecision_choleskyDecomp: A must be nonempty"
        end if

        R = A
        do k = 1, m-1
            R(k+1:m,k) = 0._qp
        end do

        do k = 1, m
            if (.NOT. (R(k,k) > 0._qp)) then
                error stop "linAlgQuadPrecision_choleskyDecomp: A is not positive definite"
            end if
            do j = k+1, m
                R(j,j:m) = R(j,j:m) - R(k,j:m)*R(k,j)/R(k,k)
            end do
            R(k,k:m) = R(k,k:m)/SQRT(R(k,k))
        end do

    end function linAlgQuadPrecision_choleskyDecomp

    function linAlgQuadPrecision_forwardSub_matRHS(L,b) result(y)
        ! Solves an mxm lower-triangular system Ly=b for y.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ly=b.
        ! The result y has the same shape as b.
        implicit none
        real(qp), intent(in) :: L(:,:), b(:,:)
        real(qp) :: y(SIZE(b,1),SIZE(b,2))
        real(qp) :: ss(SIZE(b,2))
        integer :: i,j, m,n

        m = SIZE(L,1)
        n = SIZE(b,2)

        if (SIZE(L,2) /= m) then
            error stop "linAlgQuadPrecision_forwardSub: L must be square"
        end if
        if (SIZE(b,1) /= m) then
            error stop "linAlgQuadPrecision_forwardSub: SIZE(b,1) must be equal to SIZE(L,1)"
        end if
        
        y(1,:) = b(1,:)/L(1,1)
        do i=2,m
            ss=0._qp
            do j=1,i-1
                ss = ss + L(i,j)*y(j,:)
            end do
            y(i,:) = (b(i,:) - ss)/L(i,i)
        end do
    end function linAlgQuadPrecision_forwardSub_matRHS

    function linAlgQuadPrecision_forwardSub_vecRHS(L, b) result(y)
        ! Solves an mxm lower-triangular system Ly=b for y.
        ! b is an mx1 column vector
        ! The result y has the same shape as b.
        implicit none
        real(qp), intent(in) :: L(:,:), b(:)
        real(qp) :: y(SIZE(b))
        real(qp) :: ss
        integer :: i,j, m

        m = SIZE(L,1)

        if (SIZE(L,2) /= m) then
            error stop "linAlgQuadPrecision_forwardSub: L must be square"
        end if
        if (SIZE(b) /= m) then
            error stop "linAlgQuadPrecision_forwardSub: SIZE(b) must be equal to SIZE(L,1)"
        end if
        
        y(1) = b(1)/L(1,1)
        do i=2,m
            ss=0._qp
            do j=1,i-1
                ss = ss + L(i,j)*y(j)
            end do
            y(i) = (b(i) - ss)/L(i,i)
        end do
    end function linAlgQuadPrecision_forwardSub_vecRHS

    function linAlgQuadPrecision_backSub_matRHS(U, y) result(x)
        ! Solves an mxm upper-triangular system Ux=y for x
        ! y is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ux=y.
        ! The result x has the same shape as y.
        implicit none
        real(qp), intent(in) :: U(:,:), y(:,:)
        real(qp) :: x(SIZE(y,1),SIZE(y,2))
        real(qp) :: ss(SIZE(y,2))
        integer :: i,j, m,n

        m = SIZE(U,1)
        n = SIZE(y,2)

        if (SIZE(U,2) /= m) then
            error stop "linAlgQuadPrecision_backSub: U must be square"
        end if
        if (SIZE(y,1) /= m) then
            error stop "linAlgQuadPrecision_backSub: SIZE(y,1) must be equal to SIZE(U,1)"
        end if
        
        x(m,:) = y(m,:)/U(m,m)
        do i=m-1,1,-1
            ss=0._qp
            do j=i+1,m
                ss = ss + U(i,j)*x(j,:)
            end do
            x(i,:) = (y(i,:) - ss)/U(i,i)
        end do
    end function linAlgQuadPrecision_backSub_matRHS

    function linAlgQuadPrecision_backSub_vecRHS(U, y) result(x)
        ! Solves an mxm upper-triangular system Ux=y for x
        ! y is an mx1 column vector
        ! The result x has the same shape as y.
        implicit none
        real(qp), intent(in) :: U(:,:), y(:)
        real(qp) :: x(SIZE(y))
        real(qp) :: ss
        integer :: i,j, m
        
        m = SIZE(U,1)

        if (SIZE(U,2) /= m) then
            error stop "linAlgQuadPrecision_backSub: U must be square"
        end if
        if (SIZE(y) /= m) then
            error stop "linAlgQuadPrecision_backSub: SIZE(y) must be equal to SIZE(U,1)"
        end if
        
        x(m) = y(m)/U(m,m)
        do i=m-1,1,-1
            ss=0._qp
            do j=i+1,m
                ss = ss + U(i,j)*x(j)
            end do
            x(i) = (y(i) - ss)/U(i,i)
        end do
    end function linAlgQuadPrecision_backSub_vecRHS

    function linAlgQuadPrecision_solve_matRHS(A, b) result(x)
        ! Solves the mxm system Ax=b for x.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ax=b.
        ! The result x has the same shape as b.
        implicit none
        real(qp), intent(in)  :: A(:,:), b(:,:)
        real(qp) :: x(SIZE(b,1),SIZE(b,2))
        real(qp) :: L(SIZE(A,1),SIZE(A,1)), U(SIZE(A,1),SIZE(A,1))
        integer :: p(SIZE(b,1))
        real(qp) :: Pb(SIZE(b,1),SIZE(b,2))
        integer :: i, m,n

        m = SIZE(A,1)
        n = SIZE(b,2)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_solve: A must be square"
        end if
        if (SIZE(b,1) /= m) then
            error stop "linAlgQuadPrecision_solve: SIZE(b,1) must be equal to SIZE(A,1)"
        end if

        call linAlgQuadPrecision_LU(A, L, U, p)

        ! apply permutation vector to b
        do i=1,m
            Pb(i,:) = b(p(i), :)
        end do

        ! now we are set up to solve LUx = PAx = Pb for x.
        x = linAlgQuadPrecision_forwardSub(L, Pb) ! first solve Lx' = Pb where x'=Ux
        x = linAlgQuadPrecision_backSub(U, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solve_matRHS

    function linAlgQuadPrecision_solve_vecRHS(A, b) result(x)
        ! Solves the mxm system Ax=b for x.
        ! b is an mx1 column vector
        ! The result x has the same shape as b.
        implicit none
        real(qp), intent(in)  :: A(:,:), b(:)
        real(qp) :: x(SIZE(b))
        real(qp) :: L(SIZE(A,1),SIZE(A,1)), U(SIZE(A,1),SIZE(A,1))
        integer :: p(SIZE(b))
        real(qp) :: Pb(SIZE(b))
        integer :: i, m

        m = SIZE(A,1)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_solve: A must be square"
        end if
        if (SIZE(b) /= m) then
            error stop "linAlgQuadPrecision_solve: SIZE(b) must be equal to SIZE(A,1)"
        end if

        call linAlgQuadPrecision_LU(A, L, U, p)

        ! apply permuation vector to b
        do i=1,m
            Pb(i) = b(p(i))
        end do

        ! now we are set up to solve LUx = PAx = Pb for x.
        x = linAlgQuadPrecision_forwardSub(L, Pb) ! first solve Lx' = Pb where x'=Ux
        x = linAlgQuadPrecision_backSub(U, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solve_vecRHS

    function linAlgQuadPrecision_solveSPD_matRHS(A, b) result(x)
        ! Solves the mxm system Ax=b for x when A is symmetric positive definite.
        ! b is an mxn matrix which corresponds to simultaneously solving
        ! n systems of Ax=b.
        ! The result x has the same shape as b.
        implicit none
        real(qp), intent(in) :: A(:,:)
        real(qp), intent(in) :: b(:,:)
        real(qp) :: x(SIZE(b,1),SIZE(b,2))
        real(qp) :: R(SIZE(A,1),SIZE(A,1)), Rt(SIZE(A,1),SIZE(A,1))
        integer :: m,n

        m = SIZE(A,1)
        n = SIZE(b,2)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_solveSPD: A must be square"
        end if
        if (SIZE(b,1) /= m) then
            error stop "linAlgQuadPrecision_solveSPD: SIZE(b,1) must be equal to SIZE(A,1)"
        end if

        R = linAlgQuadPrecision_choleskyDecomp(A)
        Rt = TRANSPOSE(R)

        ! now we are set up to solve LUx = Ax = b for x where L=R^T and U=R
        x = linAlgQuadPrecision_forwardSub(Rt, b) ! first solve Lx' = b where x'=Ux
        x = linAlgQuadPrecision_backSub(R, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solveSPD_matRHS

    function linAlgQuadPrecision_solveSPD_vecRHS(A, b) result(x)
        ! Solves the mxm system Ax=b for x when A is symmetric positive definite.
        ! b is an mx1 column vector
        ! The result x has the same shape as b.
        implicit none
        real(qp), intent(in) :: A(:,:)
        real(qp), intent(in) :: b(:)
        real(qp) :: x(SIZE(b))
        real(qp) :: R(SIZE(A,1),SIZE(A,1)), Rt(SIZE(A,1),SIZE(A,1))
        integer :: m

        m = SIZE(A,1)

        if (SIZE(A,2) /= m) then
            error stop "linAlgQuadPrecision_solveSPD: A must be square"
        end if
        if (SIZE(b) /= m) then
            error stop "linAlgQuadPrecision_solveSPD: SIZE(b) must be equal to SIZE(A,1)"
        end if

        R = linAlgQuadPrecision_choleskyDecomp(A)
        Rt = TRANSPOSE(R)

        ! now we are set up to solve LUx = Ax = b for x where L=R^T and U=R
        x = linAlgQuadPrecision_forwardSub(Rt, b) ! first solve Lx' = b where x'=Ux
        x = linAlgQuadPrecision_backSub(R, x) ! now solve Ux = x' for x
        
    end function linAlgQuadPrecision_solveSPD_vecRHS

end module linAlgQuadPrecision
