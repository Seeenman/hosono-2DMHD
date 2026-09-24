program test_reconstruct_faceValsWithGP
    ! program:      test_reconstruct_faceValsWithGP

    use definitions, only: max_string_length, xdir, ydir
    use assert, only: assert_close, assert_summary
    use GP, only: GP_init, GP_finalize, GP_nQuadratureMax, GP_maxRadius
    use grid, only: grid_init, grid_finalize, grid_block
    use gridBlock, only: gridBlock_t
    use initialCondition, only: initialCondition_set
    use reconstruct, only: reconstruct_faceValsWithGP
    use simulation, only: simulation_init, simulation_finalize

    implicit none

    character(len=max_string_length) :: paramfile
    integer :: exitStat
    integer :: radius
    integer :: i, j
    integer :: vv

    paramfile = "par/test_reconstruct_faceValsWithGP.par"

    call GP_init(paramfile)
    call grid_init(paramfile, GP_maxRadius, GP_nQuadratureMax)
    call simulation_init(paramfile)

    call initialCondition_set(paramfile, grid_block)

    radius = 4

    call reconstruct_faceValsWithGP(grid_block, radius)
    i = 10
    j = 10
    vv = 3
    call print_faceValGrid(grid_block, vv, i, j, radius+1)
    ! call reconstruct_faceValsWithGP(grid_block, 2)
    ! call reconstruct_faceValsWithGP(grid_block, 3)

    call simulation_finalize()
    call grid_finalize()
    call GP_finalize()

    call assert_summary(exitStat)
    call EXIT(exitStat)

contains

    subroutine print_faceValGrid(blk, vv, ic, jc, nq)
        ! purpose:   Print the cell averages and reconstructed face values of
        !            one conserved variable on the 3x3 block of cells centered
        !            on cell (ic,jc). Each cell's face values are printed just
        !            inside the face they belong to, so the two Riemann states
        !            on a shared face sit on either side of the line.
        !
        ! Inputs:    - blk (gridBlock_t) block holding U, lowerFace, upperFace
        !            - vv (integer) index of the conserved variable to print
        !            - ic, jc (integer) indices of the center cell
        !            - nq (integer) number of quadrature points per face
        ! ------------------------------------------------------------
        implicit none
        type(gridBlock_t), intent(in) :: blk
        integer, intent(in) :: vv, ic, jc, nq
        ! local variables
        integer, parameter :: vw = 12              ! width of one ES12.5 value
        integer, parameter :: labw = 16            ! room for the "U(v,i,j)=" label
        integer, parameter :: gap = 5              ! min space between side values and label
        character(len=*), parameter :: vfmt = '(ES12.5)'
        character(len=:), allocatable :: lines(:)
        character(len=32) :: label
        integer :: cw, nm, ch, width, nlines
        integer :: r, c, i, j, k, m, p
        integer :: line0, col0, tw, strt

        cw = MAX(2*vw + 2*gap + labw, nq*vw + (nq-1)*2 + 2) ! cell interior width
        nm = MAX(nq, 2)                                      ! rows for the x faces
        ch = nm + 4                                          ! cell interior height
        width = 3*(cw+1) + 1
        nlines = 3*(ch+1) + 1

        allocate(character(len=width) :: lines(nlines))
        lines = ' '

        ! grid lines
        do r = 0, 3
            lines(r*(ch+1)+1) = REPEAT('-', width)
        end do
        do r = 0, 2
            do m = 1, ch
                do c = 0, 3
                    lines(r*(ch+1)+1+m)(c*(cw+1)+1:c*(cw+1)+1) = '|'
                end do
            end do
        end do

        ! quadraturePoints(1) is the most positive, so k=1 sits at the right of
        ! the y faces and at the top of the x faces
        do r = 0, 2
            j = jc + 1 - r
            line0 = r*(ch+1) + 1
            do c = 0, 2
                i = ic - 1 + c
                col0 = c*(cw+1) + 1

                ! upper (top) and lower (bottom) y faces
                tw = nq*vw + (nq-1)*2
                strt = col0 + 1 + (cw - tw)/2
                do p = 0, nq-1
                    k = nq - p
                    write(lines(line0+1)(strt+p*(vw+2):strt+p*(vw+2)+vw-1), vfmt) blk%upperFace(vv, ydir, k, i, j)
                    write(lines(line0+ch)(strt+p*(vw+2):strt+p*(vw+2)+vw-1), vfmt) blk%lowerFace(vv, ydir, k, i, j)
                end do

                ! lower (left) and upper (right) x faces
                do k = 1, nq
                    m = line0 + 2 + (nm - nq)/2 + k
                    write(lines(m)(col0+2:col0+1+vw), vfmt) blk%lowerFace(vv, xdir, k, i, j)
                    write(lines(m)(col0+cw-vw:col0+cw-1), vfmt) blk%upperFace(vv, xdir, k, i, j)
                end do

                ! cell average in the center
                write(label, '(A,I0,A,I0,A,I0,A)') "U(", vv, ",", i, ",", j, ")="
                m = line0 + 2 + nm/2
                strt = col0 + 1 + (cw - LEN_TRIM(label))/2
                lines(m)(strt:strt+LEN_TRIM(label)-1) = TRIM(label)
                strt = col0 + 1 + (cw - vw)/2
                write(lines(m+1)(strt:strt+vw-1), vfmt) blk%U(vv, i, j)
            end do
        end do

        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') "U(", vv, ",:,:) and its face values on the 3x3 cells centered on (i,j) = (", &
            ic, ",", jc, ") with ", nq, " quadrature points per face"
        write(*,'(A)') "  i increases to the right, j increases upward"
        write(*,'(A)') "  top/bottom of each cell: upperFace/lowerFace(vv, ydir, k, i, j), k = nq..1 from left to right"
        write(*,'(A)') "  left/right of each cell: lowerFace/upperFace(vv, xdir, k, i, j), k = 1..nq from top to bottom"
        do m = 1, nlines
            write(*,'(A)') lines(m)
        end do

        deallocate(lines)
    end subroutine print_faceValGrid

end program test_reconstruct_faceValsWithGP
