module grid

    ! Module that owns the instance of the grid block custom data
    ! type for the run as well as the subroutines to initialize
    ! (allocate) and finalize (deallocate) it.
    !
    ! The grid data for a block (in a serial simulatoin the entire
    ! computational domain) is in grid_block which is of type gridBlock_t
    ! (defined in gridBlock.f90) and includes indices, domain geometry,
    ! primitive and conservative variables at each cell, face centered
    ! values, etc).

    use definitions, only: max_string_length, ndim, xdir, ydir
    use readParamFile, only: readParamFile_int, readParamFile_real
    use gridBlock, only: gridBlock_t, gridBlock_alloc, gridBlock_dealloc
    use GP, only: GP_radius, GP_nQuadrature

    implicit none

    ! The grid for this run. See gridBlock.f90 for data in this struct
    type(gridBlock_t) :: grid_block

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

        ! read values in from parameter file
        grid_block%N(xdir) = readParamFile_int(paramfile, "grid_Nx")
        grid_block%N(ydir) = readParamFile_int(paramfile, "grid_Ny")
        grid_block%domainBeg(xdir) = readParamFile_real(paramfile, "grid_xBeg")
        grid_block%domainEnd(xdir) = readParamFile_real(paramfile, "grid_xEnd")
        grid_block%domainBeg(ydir) = readParamFile_real(paramfile, "grid_yBeg")
        grid_block%domainEnd(ydir) = readParamFile_real(paramfile, "grid_yEnd")

        grid_block%NGC = GP_radius

        ! set other variables based on what was read in from the paramter file
        do i_dim=1,ndim
            grid_block%minIdx(i_dim) = 1 ! index of first guard cell
            grid_block%maxIdx(i_dim) = grid_block%N(i_dim) + 2*grid_block%NGC ! index of last guard cell
            grid_block%strtIdx(i_dim) = grid_block%minIdx(i_dim) + grid_block%NGC ! index of first real (interior) cell
            grid_block%stopIdx(i_dim) = grid_block%maxIdx(i_dim) - grid_block%NGC ! index of last real (interior) cell
            ! set dx and dy
            grid_block%dl(i_dim) = (grid_block%domainEnd(i_dim)-grid_block%domainBeg(i_dim))/grid_block%N(i_dim)
        end do

        !!!! allocate (and zero) the arrays belonging to the block, now
        !!!! that we know Nx, Ny, the guard cell count, and the
        !!!! number of quadrature points
        call gridBlock_alloc(grid_block, GP_nQuadrature)

        ! fill in grid points
        do i=grid_block%minIdx(xdir), grid_block%maxIdx(xdir)
            grid_block%x(i) = (i-grid_block%NGC-0.5)*grid_block%dl(xdir) + grid_block%domainBeg(xdir)
        end do
        do i=grid_block%minIdx(ydir), grid_block%maxIdx(ydir)
            grid_block%y(i) = (i-grid_block%NGC-0.5)*grid_block%dl(ydir) + grid_block%domainBeg(ydir)
        end do

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
        call gridBlock_dealloc(grid_block)
    end subroutine grid_finalize

end module grid
