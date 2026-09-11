module GP_CoVMatSqExp

    implicit none

contains

    subroutine GP_CovMatSqExp_(x1, x2)
        ! subroutine:   GP_CovMatSqExp_
        ! purpose:      Compute the GP covariance matrix using the squared exponential kernel
        !               k(x,y) = exp(-|x-y|^2/(2l^2)) for volume averaged quantities.
        !               The 
        ! 
        ! Inputs:       
        !               
        ! Outputs:      
        ! ------------------------------------------------------------

        implicit none
        
        real, intent(IN) :: x1(:,:), x2(:,:)

        ! local variables
        

        return
    end subroutine GP_CovMatSqExp_

end module GP_CoVMatSqExp
