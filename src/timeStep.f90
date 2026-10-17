module timeStep

    use definitions, only: xdir, ydir, nConsVars
    use gridBlock, only: gridBlock_t
    use gridState, only: gridState_t
    use convert, only: convert_cons2prim
    use boundaryConditions, only: boundaryConditions_applyBlk
    use reconstruct, only: reconstruct_faceValsWithGP
    use getFluxes, only: getFluxes_
    use GP, only: GP_maxRadius
    use detect, only: detect_ALL

    implicit none

    private

    public timeStep_advanceSolution

contains

    subroutine timeStep_advanceSolution(dt, blk, state)
        implicit none
        real, intent(in) :: dt
        type(gridBlock_t), intent(in out) :: blk
        type(gridState_t), intent(in out) :: state
        ! local variables
        integer :: i, j

        state%U = timeStep_RK3(dt, blk, state%U)
        
        !$omp parallel do collapse(2) default(none) shared(state) private(i,j)
        do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
            do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                state%V(:,i,j) = convert_cons2prim(state%U(:,i,j))
            end do
        end do
        !$omp end parallel do

        call boundaryConditions_applyBlk(blk, state%U)
        call boundaryConditions_applyBlk(blk, state%V)
        
    end subroutine timeStep_advanceSolution

    function timeStep_forwardEuler(dt, blk, U) result(Unew)
        ! purpose:   TODO
        !            
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - U (real array) the conservative variables on blk
        !              to advance, e.g., grid_state%U or the U of an RK
        !              substage
        !            
        ! Outputs:   TODO
        !            
        !            
        !            
        ! ------------------------------------------------------------
        implicit none
        real, intent(in) :: dt
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        real :: Unew(nConsVars, &
                     blk%minIdx(xdir):blk%maxIdx(xdir), &
                     blk%minIdx(ydir):blk%maxIdx(ydir))
        ! local variables
        real :: dx, dy
        integer :: i, j
        logical :: exit_mood_loop

        dx = blk%dl(xdir)
        dy = blk%dl(ydir)

        ! start each RK substage with full high order reconstruction
        ! and compute everything on the first pass of the MOOD loop
        blk%scheme = GP_maxRadius
        blk%recompute = .true.

        exit_mood_loop = .false.

        ! The MOOD loop
        do while(.not. exit_mood_loop)

            ! reconstruct face values at all quadrature points at the correct
            ! order based on blk%scheme. blk%scheme holds the radius of the
            ! GP method to be used for each face of each cell. 
            call reconstruct_faceValsWithGP(blk, U)

            ! fill in fluxes. Don't need to pass in U because this reads
            ! directly from blk%faceVals which was updated by
            ! reconstruct_faceValsWithGP using whatever substage value of U
            ! we are on.
            call getFluxes_(blk)
            ! Call FLICDdivB here

            ! only cells flagged by the MOOD loop have new fluxes. Unew
            ! keeps the previous pass's values everywhere else

            !$omp parallel do collapse(2) default(none) shared(blk,U,Unew,dt,dx,dy) private(i,j)
            do j=blk%strtIdx(ydir), blk%stopIdx(ydir)
                do i=blk%strtIdx(xdir), blk%stopIdx(xdir)
                    if (.not. blk%recompute(i,j)) CYCLE
                    ! the unsplit finite volume update
                    Unew(:,i,j) = U(:,i,j) &
                        - (dt/dx)*(blk%flux(:,xdir,i+1,j) - blk%flux(:,xdir,i,j)) &
                        - (dt/dy)*(blk%flux(:,ydir,i,j+1) - blk%flux(:,ydir,i,j))
                end do
            end do
            !$omp end parallel do

            call detect_ALL(blk, U, Unew, exit_mood_loop)

        end do

    end function timeStep_forwardEuler

    function timeStep_RK3(dt, blk, U) result(Unew)
        ! purpose:   TODO
        !            
        ! 
        ! Inputs:    - blk (gridBlock_t) the block on which to reconstruct
        !              pointwise values at cell faces
        !            - U (real array) the conservative variables on blk
        !              to advance, e.g., grid_state%U or the U of an RK
        !              substage
        !            
        ! Outputs:   TODO
        !            
        !            
        !            
        ! ------------------------------------------------------------
        implicit none
        real, intent(in) :: dt
        type(gridBlock_t), intent(in out) :: blk
        real, intent(in) :: U(nConsVars, &
                              blk%minIdx(xdir):blk%maxIdx(xdir), &
                              blk%minIdx(ydir):blk%maxIdx(ydir))
        real :: Unew(nConsVars, &
                     blk%minIdx(xdir):blk%maxIdx(xdir), &
                     blk%minIdx(ydir):blk%maxIdx(ydir))
        ! local variables
        real, dimension(nConsVars, &
            blk%minIdx(xdir):blk%maxIdx(xdir), &
            blk%minIdx(ydir):blk%maxIdx(ydir)) :: U1, U2, U3

        U1 = timeStep_forwardEuler(dt, blk, U)
        call boundaryConditions_applyBlk(blk, U1)
        U2 = timeStep_forwardEuler(dt, blk, U1)
        call boundaryConditions_applyBlk(blk, U2)
        U2 = (3*U+U2)/4.0
        call boundaryConditions_applyBlk(blk, U2)
        U3 = timeStep_forwardEuler(dt, blk, U2)
        call boundaryConditions_applyBlk(blk, U3)
        Unew = (U + 2*U3)/3.0

    end function timeStep_RK3

end module timeStep
