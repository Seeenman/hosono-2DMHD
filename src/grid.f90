module grid

    ! Module that owns the instances of the grid block and grid state
    ! custom data types for the run as well as the subroutines to
    ! initialize (allocate) and finalize (deallocate) them.
    !
    ! The grid data for a block (in a serial simulatoin the entire
    ! computational domain) is in grid_block which is of type gridBlock_t
    ! (defined in gridBlock.f90) and includes indices, domain geometry,
    ! face centered values, fluxes, etc.
    !
    ! The solution living on that block (primitive and conservative
    ! variables at each cell) is in grid_state which is of type
    ! gridState_t (defined in gridState.f90).

    use definitions, only: max_string_length, ndim, xdir, ydir
    use readParamFile, only: readParamFile_int, readParamFile_real
    use gridBlock, only: gridBlock_t, gridBlock_alloc, gridBlock_dealloc
    use gridState, only: gridState_t, gridState_alloc, gridState_dealloc

    implicit none

    ! The grid for this run. See gridBlock.f90 for data in this struct
    type(gridBlock_t) :: grid_block

    ! The solution on grid_block. See gridState.f90 for data in this struct
    type(gridState_t) :: grid_state

contains

    subroutine grid_init(paramfile, GP_maxRadius, GP_nQuadratureMax)
        ! subroutine:   grid_init
        ! Author:       Sean Riedel
        ! purpose:      To read in (from a parameter file) all grid 
        !               related parameters, and then allocate and 
        !               fill in grid related variables
        ! 
        ! Inputs:       - paramfile (character) name of file to be read
        !               - GP_maxRadius (integer) max radius of GP stencil
        !               - GP_nQuadratureMax (integer) max number of Gauss-Legendre
        !                 quadrature points per face
        !               
        ! Outputs:      - 
        ! ------------------------------------------------------------
        
        implicit none
        ! subroutine arguments
        character(len=max_string_length), intent(in) :: paramfile
        integer, intent(in) :: GP_maxRadius, GP_nQuadratureMax
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

        ! Fluxes are needed one face beyond the boundary faces of the
        ! interior (for cell centered constrained transport), so face values
        ! are reconstructed on two rings of guard cells (see reconstruct.f90).
        ! The GP stencil reaches GP_maxRadius cells along each axis, so
        ! reconstructing at strtIdx-2 needs cells down to strtIdx-2-GP_maxRadius,
        ! i.e. GP_maxRadius+2 guard cells.
        grid_block%NGC = GP_maxRadius+2

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
        call gridBlock_alloc(grid_block, GP_nQuadratureMax)
        call gridState_alloc(grid_state, grid_block)

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
        call gridState_dealloc(grid_state)
        call gridBlock_dealloc(grid_block)
        write(*,*) "=============================================================="
        write(*,*) "Grid deallocated."
        write(*,*) "=============================================================="
    end subroutine grid_finalize

end module grid
