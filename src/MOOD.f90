module MOOD
    ! module with MOOD (Multidimensional Optimal Order Detection)
    ! related subroutines.

    use definitions, only: ndim
    use gridBlock, only: gridBlock_t
    use GP, only: GP_maxRadius

    implicit none

    private

    public :: MOOD_loop

contains

    subroutine MOOD_loop(blk)
        ! The GP-MOOD loop that determines
        ! which order to use to update the solution at each cell
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in out) :: blk
        ! local variables
        integer :: 
        integer :: i,j

        ! unlimited GP reconstruction of U

    end subroutine MOOD_loop

end module MOOD
