module output

    ! Writes simulation data to disk.
    !
    ! Note: this module was adapted from a 3D code. The z-direction
    ! datasets have been dropped and the field datasets are now 2D.
    !
    ! The grid data this module writes is passed in through the argument
    ! list rather than being use-associated from the grid module, so that
    ! this module has no knowledge of where that data lives. It arrives
    ! bundled as a gridBlock_t.

    use hdf5
    use definitions, only: ndim, xdir, ydir, max_string_length, &
                           dens_var, velx_var, vely_var, velz_var, &
                           magx_var, magy_var, magz_var, &
                           pres_var, eint_var
    use gridBlock, only: gridBlock_t
    use simulation, only: sim_dataFileBaseName, sim_outputFreqStep, sim_outputFreqTime, &
                          sim_outputHdf5, sim_outputAscii, sim_gamma

    implicit none

    private

    public :: output_write

contains

    subroutine output_write(nStep, t, dt, lastOutputStep, lastOutputTime, outputCounter, forceOutput, blk)
        ! purpose:      Decide whether it is time to write an output file,
        !               and write one if so
        !
        ! Inputs:       - nStep (integer) current time step number
        !               - t/dt (real) current time and time step size
        !               - lastOutputStep/lastOutputTime (integer/real) step and
        !                 time at which the previous output was written
        !               - outputCounter (integer) number of the next output file
        !               - forceOutput (logical) write regardless of the output
        !                 frequency (used for the initial condition, for example)
        !               - blk (gridBlock_t) the block of the grid to write,
        !                 carrying its own indices, geometry, and fluid state
        !
        ! Note:         blk%V must have up to date guard cells, because the
        !               divergence of B is evaluated with a centered stencil
        !               that reaches one cell outside the interior.
        !
        ! Outputs:      - lastOutputStep/lastOutputTime/outputCounter are updated
        !                 whenever a file is written
        ! ------------------------------------------------------------
        implicit none
        real, intent(in) :: t, dt
        real, intent(inout) :: lastOutputTime
        integer, intent(in) :: nStep
        integer, intent(inout) :: lastOutputStep, outputCounter
        logical, intent(in) :: forceOutput
        type(gridBlock_t), intent(in) :: blk
        ! local variables
        character(len=max_string_length) :: outputfile
        character(len=5) :: counterChar
        logical :: writenow

        ! check to make sure at least one of sim_outputFreqStep and sim_outputFreqTime is positive
        if ((sim_outputFreqTime<=0) .and. (sim_outputFreqStep<=0)) then
            write(*,*) "=========================================================================="
            write(*,*) "No ouptputs will be written because both sim_outputFreqStep and"
            write(*,*) "sim_outputFreqTime are not positive. Please set at least one of them"
            write(*,*) "to be greater than zero."
            write(*,*) "=========================================================================="
            stop
        end if

        ! check to see if it is time to write an output file
        writenow = .false.
        if ((nStep-lastOutputStep == sim_outputFreqStep) .and. (sim_outputFreqStep>0)) then
            writenow = .true.
        else if ((t-lastOutputTime >= sim_outputFreqTime) .and. (sim_outputFreqTime>0)) then
            writenow = .true.
        end if

        if (writenow .or. forceOutput) then

            ! convert counter number to a character
            write(counterChar, '(i5.5)') outputCounter

            if (sim_outputAscii) then
                ! file name for ascii output
                ! outputfile = trim(sim_dataFileBaseName)//'_'//trim(counterChar)//'.dat'
                ! call output_writeAscii(nStep, t, outputCounter, outputfile)
                write(*,*)'Ascii output not supported right now'
                stop
            end if
            if (sim_outputHdf5) then
                ! file name for hdf5 output
                outputfile = trim(sim_dataFileBaseName)//'_'//trim(counterChar)//'.h5'
                call output_writeHdf5(nStep, t, dt, outputCounter, trim(outputfile), blk)
            end if

            ! write a message to stdout
            write(*,*)''
            write(*,*)' Output no.',outputCounter, 'has been written      '
            write(*,*)'================================================='
            write(*,*)'   Steps      Time              dt               '
            write(*,*)'================================================='
            write(*,*)''

            ! update last output time, last output step, and outputCounter
            lastOutputStep = nStep
            lastOutputTime = outputCounter*sim_outputFreqTime
            outputCounter = outputCounter + 1

        end if

    end subroutine output_write

    subroutine output_writeHdf5(nStep, t, dt, outputCounter, outputfile, blk)
        ! purpose:      Write one hdf5 output file
        !
        ! Inputs:       - outputfile (character) name of the file to write
        !               - all other arguments are as described in output_write
        !
        ! Outputs:      - none (writes a file to disk)
        ! ------------------------------------------------------------
        implicit none
        real, intent(in) :: t, dt
        integer, intent(in) :: nStep, outputCounter
        character(len=*), intent(in) :: outputfile
        type(gridBlock_t), intent(in) :: blk
        ! local variables
        ! Buffers holding the interior of each field, contiguous and
        ! guard-cell free, ready to hand to hdf5. Allocatable rather than
        ! automatic so that where they live does not depend on whether the
        ! compiler puts large local arrays on the stack or the heap.
        real, allocatable, dimension(:,:) :: dens, velx, vely, velz, pres, eint, gama, ener
        real, allocatable, dimension(:,:) :: magp, divb, magx, magy, magz
        integer :: error, space_rank
        integer :: i_dim, ii, jj, offsets(ndim)
        integer(HSIZE_T) :: data_dims_1d(1), data_dims_2d(2)
        integer(HID_T) :: file_id, dspace_id, &
                          dset_id_outputCounter, dset_id_t, dset_id_dt, dset_id_nStep, &
                          dset_id_xmin, dset_id_xmax, dset_id_nx, dset_id_velx, &
                          dset_id_ymin, dset_id_ymax, dset_id_ny, dset_id_vely, &
                          dset_id_velz, &
                          dset_id_dens, dset_id_pres, dset_id_eint, dset_id_gama, &
                          dset_id_ener, &
                          dset_id_magp, dset_id_divb, dset_id_magx, dset_id_magy, dset_id_magz

        !~~~! allocate the field buffers
        allocate(dens(blk%N(xdir), blk%N(ydir)))
        allocate(velx(blk%N(xdir), blk%N(ydir)))
        allocate(vely(blk%N(xdir), blk%N(ydir)))
        allocate(velz(blk%N(xdir), blk%N(ydir)))
        allocate(pres(blk%N(xdir), blk%N(ydir)))
        allocate(eint(blk%N(xdir), blk%N(ydir)))
        allocate(gama(blk%N(xdir), blk%N(ydir)))
        allocate(ener(blk%N(xdir), blk%N(ydir)))
        allocate(magp(blk%N(xdir), blk%N(ydir)))
        allocate(divb(blk%N(xdir), blk%N(ydir)))
        allocate(magx(blk%N(xdir), blk%N(ydir)))
        allocate(magy(blk%N(xdir), blk%N(ydir)))
        allocate(magz(blk%N(xdir), blk%N(ydir)))

        !=! open hdf5 interface
        call h5open_f(error)

        !++! open the file
        call h5fcreate_f(outputfile, H5F_ACC_TRUNC_F, file_id, error) !! H5F_ACC_TRUNC_F overwrites the file if it already exists
        if (error /= 0) then
            write(*,*) "=========================================================================="
            write(*,*) "Could not create the hdf5 output file: ", outputfile
            write(*,*) "Check that the directory it lives in exists and is writable."
            write(*,*) "(the directory comes from sim_dataFileBaseName in the parameter file)"
            write(*,*) "=========================================================================="
            stop
        end if

        !:::! open dataspace for outputCounter, t, and nStep
        space_rank = 1 ! number of dimensions in the data space
        data_dims_1d(1) = 1
        call h5screate_simple_f(space_rank, data_dims_1d, dspace_id, error)

        !\\\\! create datasets for outputCounter, t, and nStep
        call h5dcreate_f(file_id, "output_number", H5T_NATIVE_INTEGER, dspace_id, dset_id_outputCounter, error)
        call h5dcreate_f(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id_t, error)
        call h5dcreate_f(file_id, "dt", H5T_NATIVE_DOUBLE, dspace_id, dset_id_dt, error)
        call h5dcreate_f(file_id, "step_number", H5T_NATIVE_INTEGER, dspace_id, dset_id_nStep, error)

        !*****! write datasets
        call h5dwrite_f(dset_id_outputCounter, H5T_NATIVE_INTEGER, outputCounter, data_dims_1d, error)
        call h5dwrite_f(dset_id_t, H5T_NATIVE_DOUBLE, t, data_dims_1d, error)
        call h5dwrite_f(dset_id_dt, H5T_NATIVE_DOUBLE, dt, data_dims_1d, error)
        call h5dwrite_f(dset_id_nStep, H5T_NATIVE_INTEGER, nStep, data_dims_1d, error)

        !////! close datasets for outputCounter, t, and nStep
        call h5dclose_f(dset_id_outputCounter, error)
        call h5dclose_f(dset_id_t, error)
        call h5dclose_f(dset_id_dt, error)
        call h5dclose_f(dset_id_nStep, error)

        !:::! close dataspace for outputCounter, t, and nStep
        call h5sclose_f(dspace_id, error)

        !---! open dataspace for primitive variables
        space_rank = 2 ! number of dimensions in the data space
        data_dims_2d(1) = blk%N(xdir)
        data_dims_2d(2) = blk%N(ydir)
        call h5screate_simple_f(space_rank, data_dims_2d, dspace_id, error)

        !\\\\! create datasets for primitive variables
        call h5dcreate_f(file_id, "dens", H5T_NATIVE_DOUBLE, dspace_id, dset_id_dens, error)
        call h5dcreate_f(file_id, "velx", H5T_NATIVE_DOUBLE, dspace_id, dset_id_velx, error)
        call h5dcreate_f(file_id, "vely", H5T_NATIVE_DOUBLE, dspace_id, dset_id_vely, error)
        call h5dcreate_f(file_id, "velz", H5T_NATIVE_DOUBLE, dspace_id, dset_id_velz, error)
        call h5dcreate_f(file_id, "magx", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magx, error)
        call h5dcreate_f(file_id, "magy", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magy, error)
        call h5dcreate_f(file_id, "magz", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magz, error)
        call h5dcreate_f(file_id, "pres", H5T_NATIVE_DOUBLE, dspace_id, dset_id_pres, error)
        call h5dcreate_f(file_id, "eint", H5T_NATIVE_DOUBLE, dspace_id, dset_id_eint, error)
        call h5dcreate_f(file_id, "gama", H5T_NATIVE_DOUBLE, dspace_id, dset_id_gama, error)
        call h5dcreate_f(file_id, "ener", H5T_NATIVE_DOUBLE, dspace_id, dset_id_ener, error)
        call h5dcreate_f(file_id, "magp", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magp, error)
        call h5dcreate_f(file_id, "divb", H5T_NATIVE_DOUBLE, dspace_id, dset_id_divb, error)

        ! This gets ride of the warning:
        ! "Fortran runtime warning: An array temporary was created"
        ! because we are making the copies of the arrays explicitly
        ! before passing them to the hdf5 subroutine.
        ! However, this does not remove the undesired (and presumabley slow)
        ! behavior that the warning is warning us about.

        dens = blk%V(dens_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        velx = blk%V(velx_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        vely = blk%V(vely_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        velz = blk%V(velz_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        magx = blk%V(magx_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        magy = blk%V(magy_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        magz = blk%V(magz_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        pres = blk%V(pres_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        eint = blk%V(eint_var, blk%strtIdx(xdir):blk%stopIdx(xdir), blk%strtIdx(ydir):blk%stopIdx(ydir))
        ! this code carries a single adiabatic index rather than a per-cell one
        gama = sim_gamma

        ! derived quantities
        ener = dens*(velx**2+vely**2+velz**2)/2 + dens*eint ! total energy (hydro)

        ! derived MHD quantities
        magp = (magx**2 + magy**2 + magz**2)/2 ! magnetic pressure
        ener = ener + magp ! total energy (MHD)
        divb = 0.0 ! divergence of B
        do i_dim=1,ndim
            offsets = 0
            offsets(i_dim) = 1
            ii = offsets(xdir)
            jj = offsets(ydir)
            divb = divb + 1/(2*blk%dl(i_dim))*&
                    (blk%V(magx_var+i_dim-1,&
                           blk%strtIdx(xdir)+ii:blk%stopIdx(xdir)+ii,&
                           blk%strtIdx(ydir)+jj:blk%stopIdx(ydir)+jj)&
                   - blk%V(magx_var+i_dim-1,&
                           blk%strtIdx(xdir)-ii:blk%stopIdx(xdir)-ii,&
                           blk%strtIdx(ydir)-jj:blk%stopIdx(ydir)-jj))
        end do

        !*****! write data to datasets
        call h5dwrite_f(dset_id_dens, H5T_NATIVE_DOUBLE, dens, data_dims_2d, error)
        call h5dwrite_f(dset_id_velx, H5T_NATIVE_DOUBLE, velx, data_dims_2d, error)
        call h5dwrite_f(dset_id_vely, H5T_NATIVE_DOUBLE, vely, data_dims_2d, error)
        call h5dwrite_f(dset_id_velz, H5T_NATIVE_DOUBLE, velz, data_dims_2d, error)
        call h5dwrite_f(dset_id_magx, H5T_NATIVE_DOUBLE, magx, data_dims_2d, error)
        call h5dwrite_f(dset_id_magy, H5T_NATIVE_DOUBLE, magy, data_dims_2d, error)
        call h5dwrite_f(dset_id_magz, H5T_NATIVE_DOUBLE, magz, data_dims_2d, error)
        call h5dwrite_f(dset_id_pres, H5T_NATIVE_DOUBLE, pres, data_dims_2d, error)
        call h5dwrite_f(dset_id_eint, H5T_NATIVE_DOUBLE, eint, data_dims_2d, error)
        call h5dwrite_f(dset_id_gama, H5T_NATIVE_DOUBLE, gama, data_dims_2d, error)
        call h5dwrite_f(dset_id_ener, H5T_NATIVE_DOUBLE, ener, data_dims_2d, error)
        call h5dwrite_f(dset_id_magp, H5T_NATIVE_DOUBLE, magp, data_dims_2d, error)
        call h5dwrite_f(dset_id_divb, H5T_NATIVE_DOUBLE, divb, data_dims_2d, error)

        !////! close datasets for primitive variables
        call h5dclose_f(dset_id_dens, error)
        call h5dclose_f(dset_id_velx, error)
        call h5dclose_f(dset_id_vely, error)
        call h5dclose_f(dset_id_velz, error)
        call h5dclose_f(dset_id_magx, error)
        call h5dclose_f(dset_id_magy, error)
        call h5dclose_f(dset_id_magz, error)
        call h5dclose_f(dset_id_pres, error)
        call h5dclose_f(dset_id_eint, error)
        call h5dclose_f(dset_id_gama, error)
        call h5dclose_f(dset_id_ener, error)
        call h5dclose_f(dset_id_magp, error)
        call h5dclose_f(dset_id_divb, error)

        !---! close dataspace for primitive variables
        call h5sclose_f(dspace_id, error)

        !:::! open dataspace for grid parameters
        space_rank = 1 ! number of dimensions in the data space
        data_dims_1d(1) = 1
        call h5screate_simple_f(space_rank, data_dims_1d, dspace_id, error)

        !\\\\! create datasets for xmin, xmax, ymin, ymax, nx, ny
        call h5dcreate_f(file_id, "xmin", H5T_NATIVE_DOUBLE, dspace_id, dset_id_xmin, error)
        call h5dcreate_f(file_id, "xmax", H5T_NATIVE_DOUBLE, dspace_id, dset_id_xmax, error)
        call h5dcreate_f(file_id, "ymin", H5T_NATIVE_DOUBLE, dspace_id, dset_id_ymin, error)
        call h5dcreate_f(file_id, "ymax", H5T_NATIVE_DOUBLE, dspace_id, dset_id_ymax, error)
        call h5dcreate_f(file_id, "nx", H5T_NATIVE_INTEGER, dspace_id, dset_id_nx, error)
        call h5dcreate_f(file_id, "ny", H5T_NATIVE_INTEGER, dspace_id, dset_id_ny, error)

        !*****! write datasets
        call h5dwrite_f(dset_id_xmin, H5T_NATIVE_DOUBLE, blk%domainBeg(xdir), data_dims_1d, error)
        call h5dwrite_f(dset_id_xmax, H5T_NATIVE_DOUBLE, blk%domainEnd(xdir), data_dims_1d, error)
        call h5dwrite_f(dset_id_ymin, H5T_NATIVE_DOUBLE, blk%domainBeg(ydir), data_dims_1d, error)
        call h5dwrite_f(dset_id_ymax, H5T_NATIVE_DOUBLE, blk%domainEnd(ydir), data_dims_1d, error)
        call h5dwrite_f(dset_id_nx, H5T_NATIVE_INTEGER, blk%N(xdir), data_dims_1d, error)
        call h5dwrite_f(dset_id_ny, H5T_NATIVE_INTEGER, blk%N(ydir), data_dims_1d, error)

        !////! close datasets for xmin, xmax, ymin, ymax, nx, ny
        call h5dclose_f(dset_id_xmin, error)
        call h5dclose_f(dset_id_xmax, error)
        call h5dclose_f(dset_id_ymin, error)
        call h5dclose_f(dset_id_ymax, error)
        call h5dclose_f(dset_id_nx, error)
        call h5dclose_f(dset_id_ny, error)

        !:::! close dataspace for grid parameters
        call h5sclose_f(dspace_id, error)

        !***! flush all data to disk before closing (helps prevent corruption on interrupt)
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        !++! close the file
        call h5fclose_f(file_id, error)

        !=! close hdf5 interface
        call h5close_f(error)

        !~~~! deallocate the field buffers
        deallocate(dens)
        deallocate(velx)
        deallocate(vely)
        deallocate(velz)
        deallocate(pres)
        deallocate(eint)
        deallocate(gama)
        deallocate(ener)
        deallocate(magp)
        deallocate(divb)
        deallocate(magx)
        deallocate(magy)
        deallocate(magz)

    end subroutine output_writeHdf5

end module output
