module choleskyDecomp

    implicit none

contains

    function choleskyDecomp_(A) result(R)
        ! cholesky decomposition of a real symmetric
        ! positive definite matrix
        !
        ! returns R in the upper triangle of A
        ! where the cholesky decomposition of A is
        ! A = R'R 
        implicit none
        real, intent(in)  :: A(:,:)
        real :: R(SIZE(A,1),SIZE(A,1))
        integer :: k,j,m
        
        m = SIZE(A,1)
        R = A
        do k=1,m
            do j=k+1,m
                R(j,j:m) = R(j,j:m) - R(k,j:m)*R(k,j)/R(k,k)
            end do
            R(k,k:m) = R(k,k:m)/SQRT(R(k,k))
        end do

    end function choleskyDecomp_

end module choleskyDecomp
