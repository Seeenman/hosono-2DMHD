program test_reconstruct_faceValsWithGP
    ! program:      test_reconstruct_faceValsWithGP

    use definitions, only: max_string_length, xdir, ydir
    use assert, only: assert_close, assert_summary
    use GP, only: GP_init, GP_finalize, GP_nQuadratureMax, GP_maxRadius
    use grid, only: grid_init, grid_finalize, grid_block
    use initialCondition, only: initialCondition_set
    use reconstruct, only: reconstruct_faceValsWithGP
    use simulation, only: simulation_init, simulation_finalize

    implicit none

    character(len=max_string_length) :: paramfile
    integer :: exitStat
    integer :: radius
    integer :: k, i, j
    integer :: vv, dir

    paramfile = "par/test_reconstruct_faceValsWithGP.par"

    call GP_init(paramfile)
    call grid_init(paramfile, GP_maxRadius, GP_nQuadratureMax)
    call simulation_init(paramfile)

    call initialCondition_set(paramfile, grid_block)

    radius = 1

    call reconstruct_faceValsWithGP(grid_block, radius)
    i = 10
    j = 10
    vv = 3
    write(*,'(A,I0,A,I0,A,I0,A,ES13.5)') "U(", vv, "," , i, ",", j, ") = ", grid_block%U(vv,i,j)
    write(*,'(A,I0,A,I0,A,I0,A,ES13.5)') "U(", vv, "," , i+1, ",", j, ") = ", grid_block%U(vv,i+1,j)
    write(*,'(A,I0,A,I0,A,I0,A,ES13.5)') "U(", vv, "," , i-1, ",", j, ") = ", grid_block%U(vv,i-1,j)
    write(*,'(A,I0,A,I0,A,I0,A,ES13.5)') "U(", vv, "," , i, ",", j+1, ") = ", grid_block%U(vv,i,j+1)
    write(*,'(A,I0,A,I0,A,I0,A,ES13.5)') "U(", vv, "," , i, ",", j-1, ") = ", grid_block%U(vv,i,j-1)
    write(*,*) "Face(var,dir,quadPoint,i,j)"
    do k=1,radius+1
        do dir=1,2
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES13.5)') "lowerFace(", vv, ", ", dir, ", ", k, ", ", i, ", ", j, ") = ", grid_block%lowerFace(vv, dir, k, i,j)
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES13.5)') "upperFace(", vv, ", ", dir, ", ", k, ", ", i, ", ", j, ") = ", grid_block%upperFace(vv, dir, k, i,j)
        end do
    end do
    ! call reconstruct_faceValsWithGP(grid_block, 2)
    ! call reconstruct_faceValsWithGP(grid_block, 3)

    call simulation_finalize()
    call grid_finalize()
    call GP_finalize()

    call assert_summary(exitStat)
    call EXIT(exitStat)

end program test_reconstruct_faceValsWithGP
