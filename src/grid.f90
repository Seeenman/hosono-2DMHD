module grid

    ! grid data and subroutines to initialize and finalize
    ! global grid arrays

    use definitions, only: max_string_length, ndim, xdir, ydir, nConsVars, nPrimVars
    use readParamFile, only: readParamFile_int, readParamFile_real

    implicit none

    ! number of guard cells padding the domain in each direction
    integer, allocatable :: grid_NGC

    ! Gaussian Process Stencil Radius (GPR)
    integer :: grid_GPR

    ! Quadrature
    integer :: grid_nquad ! number of face quadrature points
    real, allocatable :: grid_quadPoints(:), grid_quadWeights(:)

    ! The grid itself
    integer, allocatable :: grid_N(:), grid_strtIdx(:), grid_stopIdx(:), &
                                  grid_minIdx(:), grid_maxIdx(:)
    real, allocatable :: grid_beg(:), grid_end(:), grid_dl(:)
    real, allocatable :: grid_x(:) ! the values of the grid in the x direction
    real, allocatable :: grid_y(:) ! the values of the grid in the y direction

    ! conservative and primitive variables
    ! dimensions correspond to (variable, xcoordinate, ycoordinate)
    real, allocatable :: grid_U(:,:,:), grid_V(:,:,:)

    ! Pointwise conservaitve variables reconstructed at face quadrature points
    ! (variable, direction, quadrature_point, xcoordinate, ycoordinate)
    real, allocatable :: grid_lowerFace(:,:,:,:,:) 
    real, allocatable :: grid_upperFace(:,:,:,:,:)
    ! Note that these are indexed in terms of the cell for which they were
    ! calculated. For example grid_lowerFace(:,xdir,1,i,j)
    ! is the right Riemann state for the x-direction
    ! Riemann problem at the first gaussian quadrature point
    ! at interface (i-1/2,j). Similarly, grid_lowerFace(:,ydir,2,i,j)
    ! is the left Riemann state for the y-direction
    ! Riemann problem at the second gaussian quadrature point
    ! at interface (i,j+1/2)

    ! Flux values from local Riemann problems at face quadrature points
    ! (variable, direction, quadrature_point, xcoordinate, ycoordinate)
    real, allocatable :: grid_flux(:,:,:,:,:)
    ! Note that, for example, 
    ! grid_flux(:,xdir,2,i,j  )=F_{i-1/2, j    } at the  first gaussian quadrature point
    ! grid_flux(:,ydir,1,i,j  )=G_{i,     j-1/2} at the second gaussian quadrature point
    ! grid_flux(:,xdir,1,i+1,j)=F_{i+1/2, j    } at the first gaussian quadrature point
    ! et cetera

contains

    subroutine grid_init(paramfile)
        ! subroutine:   grid_init
        ! Author:       Sean Riedel
        ! purpose:      To read in (from a parameter file) all grid 
        !               related parameters, and then allocate and 
        !               fill in grid related variables
        ! 
        ! Inputs:       - paramfile (character) name of file to be read
        !               
        ! Outputs:      - none
        ! ------------------------------------------------------------
        
        implicit none
        ! subroutine arguments
        character(len=max_string_length), intent(in) :: paramfile
        ! local variables
        integer :: i, i_dim

        write(*,*) "=============================================================="
        write(*,*) "Initializing grid"
        write(*,*) "--------------------------------------------------------------"
        
        ! allocate arrays to be compatible with up to a 3D simulation
        allocate(grid_N(ndim))
        allocate(grid_minIdx(ndim)) ! index of first guard cell
        allocate(grid_maxIdx(ndim)) ! index of last guard cell
        allocate(grid_strtIdx(ndim)) ! index of first interior cell (first non-guard-cell)
        allocate(grid_stopIdx(ndim)) ! index of last interior cell (first non-guard-cell)
        allocate(grid_beg(ndim)) ! lower bounds of computational domain
        allocate(grid_end(ndim)) ! upper bounds of computatoinal domain
        allocate(grid_dl(ndim)) ! will hold dx and dy
        
        ! read values in from parameter file
        grid_N(xdir) = readParamFile_int(paramfile, "grid_Nx") 
        grid_N(ydir) = readParamFile_int(paramfile, "grid_Ny") 
        grid_GPR = readParamFile_int(paramfile, "grid_GPR") 
        grid_beg(xdir) = readParamFile_real(paramfile, "grid_xBeg")
        grid_end(xdir) = readParamFile_real(paramfile, "grid_xEnd")
        grid_beg(ydir) = readParamFile_real(paramfile, "grid_yBeg")
        grid_end(ydir) = readParamFile_real(paramfile, "grid_yEnd")

        grid_NGC = grid_GPR

        ! set other variables based on what was read in from the paramter file
        do i_dim=1,ndim
            grid_minIdx(i_dim) = 1 ! index of first guard cell
            grid_maxIdx(i_dim) = grid_N(i_dim) + 2*grid_NGC ! index of last guard cell
            grid_strtIdx(i_dim) = grid_minIdx(i_dim) + grid_NGC ! index of first real (interior) cell
            grid_stopIdx(i_dim) = grid_maxIdx(i_dim) - grid_NGC ! index of last real (interior) cell
            grid_dl(i_dim) = (grid_end(i_dim)-grid_beg(i_dim))/grid_N(i_dim) ! set dx and dy
        end do

        ! allocate space for grids now that we have Nx and Ny
        allocate(grid_x(grid_minIdx(xdir):grid_maxIdx(xdir)))
        grid_x = 0.0 ! zero out 
        allocate(grid_y(grid_minIdx(ydir):grid_maxIdx(ydir)))
        grid_y = 0.0 ! zero out 
        
        ! fill in grid points
        do i=grid_minIdx(xdir), grid_maxIdx(xdir)
            grid_x(i) = (i-grid_NGC-0.5)*grid_dl(xdir) + grid_beg(xdir)
        end do
        do i=grid_minIdx(ydir), grid_maxIdx(ydir)
            grid_y(i) = (i-grid_NGC-0.5)*grid_dl(ydir) + grid_beg(ydir)
        end do

        !!!! Set quadrature coordinates and quadrature weights
        !!!! based on grid_GPR
        grid_nquad = grid_GPR+1
        allocate(grid_quadPoints(grid_nquad))
        allocate(grid_quadWeights(grid_nquad))
        if (grid_nquad == 2) then
            ! GP spatial order of accuracy = 2*1+1 = 3
            ! use 4th order, 2 point quadrature rule
            grid_quadPoints(1) = 1.0/2.0/SQRT(3.0)
            grid_quadPoints(2) = -grid_quadPoints(1)
            grid_quadWeights(1) = 1.0/2.0
            grid_quadWeights(2) = 1.0/2.0
        else if (grid_nquad == 3) then
            ! GP spatial order of accuracy = 2*2+1 = 5
            ! use 6th order, 3 point quadrature rule
            grid_quadPoints(1) = 1.0/2.0*SQRT(3.0/5.0)
            grid_quadPoints(2) = 0.0
            grid_quadPoints(3) = -grid_quadPoints(1)
            grid_quadWeights(1) = 5.0/18.0
            grid_quadWeights(2) = 8.0/18.0
            grid_quadWeights(3) = 5.0/18.0
        else if (grid_nquad == 4) then
            ! GP spatial order of accuracy = 2*3+1 = 7
            ! use 8th order, 4 point quadrature rule
            grid_quadPoints(1) = 1.0/2.0*SQRT(3.0/7.0+2.0/7.0*SQRT(6.0/5.0))
            grid_quadPoints(2) = 1.0/2.0*SQRT(3.0/7.0-2.0/7.0*SQRT(6.0/5.0))
            grid_quadPoints(3) = -grid_quadPoints(2)
            grid_quadPoints(4) = -grid_quadPoints(1)
            grid_quadWeights(1) = (18.0-SQRT(30.0))/72.0
            grid_quadWeights(2) = (18.0+SQRT(30.0))/72.0
            grid_quadWeights(3) = (18.0+SQRT(30.0))/72.0
            grid_quadWeights(4) = (18.0-SQRT(30.0))/72.0
        end if
        
        !!!! allocate other grid data !!!!
        !! conservative variables
        allocate(grid_U(nConsVars,& 
                        grid_minIdx(xdir):grid_maxIdx(xdir),&
                        grid_minIdx(ydir):grid_maxIdx(ydir)))
        !! primitive variables
        allocate(grid_V(nPrimVars,& 
                        grid_minIdx(xdir):grid_maxIdx(xdir),&
                        grid_minIdx(ydir):grid_maxIdx(ydir)))
        !! lower face reconstructions
        allocate(grid_lowerFace(nConsVars, ndim, grid_nquad, &
                                grid_minIdx(xdir):grid_maxIdx(xdir),&
                                grid_minIdx(ydir):grid_maxIdx(ydir)))
        !! upper face reconstructions
        allocate(grid_upperFace(nConsVars, ndim, grid_nquad, &
                                grid_minIdx(xdir):grid_maxIdx(xdir),&
                                grid_minIdx(ydir):grid_maxIdx(ydir)))
        !! Fluxes
        allocate(grid_flux(nConsVars, ndim, grid_nquad, &
                           grid_minIdx(xdir):grid_maxIdx(xdir),&
                           grid_minIdx(ydir):grid_maxIdx(ydir)))

        !!! zero it all of these out !!!
        grid_U = 0.0 
        grid_V = 0.0 
        grid_lowerFace = 0.0 
        grid_upperFace = 0.0 
        grid_flux = 0.0

        write(*,*) "--------------------------------------------------------------"
        write(*,*) "Grid initialized"
        write(*,*) "=============================================================="
        
    end subroutine grid_init

    subroutine grid_finalize()
        ! subroutine:   grid_finalize
        ! Author:       Sean Riedel
        ! purpose:      To deallocate all grid related variables
        ! 
        ! Inputs:       - none
        !               
        ! Outputs:      - none
        ! ------------------------------------------------------------
        implicit none
        deallocate(grid_N)
        deallocate(grid_strtIdx)
        deallocate(grid_stopIdx)
        deallocate(grid_minIdx)
        deallocate(grid_maxIdx)
        deallocate(grid_beg)
        deallocate(grid_end)
        deallocate(grid_dl)
        deallocate(grid_U)
        deallocate(grid_V)
        deallocate(grid_lowerFace)
        deallocate(grid_upperFace)
        deallocate(grid_flux)
        deallocate(grid_x)
        deallocate(grid_y)
        deallocate(grid_quadPoints)
        deallocate(grid_quadWeights)
    end subroutine grid_finalize

end module grid

