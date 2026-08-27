module output

    use hdf5
    use definitions, only: ndim
    use simulation, only: sim_dataFileBaseName, sim_outputFreqStep, sim_outputFreqTime, &
                        sim_outputHdf5, sim_outputAscii
    use grid, only: grid_minIdx, grid_maxIdx, &
                    grid_strtIdx, grid_stopIdx, &
                    grid_beg, grid_end, &
                    grid_NGC, grid_N, &
                    grid_x, grid_y, grid_z, &
                    grid_V, grid_dl

    implicit none

    private

    public :: output_write

contains

    subroutine output_write(nStep, t, dt, lastOutputStep, lastOutputTime, outputCounter, forceOutput)
        implicit none
        real, intent(in) :: t, dt
        real, intent(inout) :: lastOutputTime
        integer, intent(in) :: nStep
        integer, intent(inout) :: lastOutputStep, outputCounter
        logical, intent(in) :: forceOutput
        character(len=MAX_STRING_LEN) :: outputfile, counterChar
        ! character(len=5) :: counterChar
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
            write(counterChar, '(i5)') outputCounter + 10000

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
                call output_writeHdf5(nStep, t, dt, outputCounter, outputfile)
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

    subroutine output_writeHdf5(nStep, t, dt, outputCounter, outputfile)
        implicit none
        real, intent(in) :: t, dt
        integer, intent(in) :: nStep, outputCounter
        character(len=*), intent(in) :: outputfile
        ! local variables
        real, dimension(grid_N(XDIR), grid_N(YDIR), grid_N(ZDIR)) :: dens, velx, vely, velz, pres, eint, gama, ener
        integer :: error, space_rank
        integer(HSIZE_T) :: data_dims_1d(1), data_dims_3d(3)
        integer(HID_T) :: file_id, dspace_id, &
                          dset_id_outputCounter, dset_id_t, dset_id_dt, dset_id_nStep, &
                          dset_id_xmin, dset_id_xmax, dset_id_nx, dset_id_velx, &
                          dset_id_ymin, dset_id_ymax, dset_id_ny, dset_id_vely, &
                          dset_id_zmin, dset_id_zmax, dset_id_nz, dset_id_velz, &
                          dset_id_dens, dset_id_pres, dset_id_eint, dset_id_gama, &
                          dset_id_ener
#ifdef MHD
        real, dimension(grid_N(XDIR), grid_N(YDIR), grid_N(ZDIR)) :: magp, divb, magx, magy, magz
        integer(HID_T) :: dset_id_magp, dset_id_divb, dset_id_magx, dset_id_magy, dset_id_magz
        integer :: i_dim, ii, jj, kk, offsets(NDIM)
#endif

        !=! open hdf5 interface
        call h5open_f(error)

        !++! open the file
        call h5fcreate_f(outputfile, H5F_ACC_TRUNC_F, file_id, error) !! H5F_ACC_TRUNC_F overwrites the file if it already exists

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
        space_rank = 3 ! number of dimensions in the data space
        data_dims_3d(1) = grid_N(XDIR)
        data_dims_3d(2) = grid_N(YDIR)
        data_dims_3d(3) = grid_N(ZDIR)
        call h5screate_simple_f(space_rank, data_dims_3d, dspace_id, error)

        !\\\\! create datasets for primitive variables
        call h5dcreate_f(file_id, "dens", H5T_NATIVE_DOUBLE, dspace_id, dset_id_dens, error)
        call h5dcreate_f(file_id, "velx", H5T_NATIVE_DOUBLE, dspace_id, dset_id_velx, error)
        call h5dcreate_f(file_id, "vely", H5T_NATIVE_DOUBLE, dspace_id, dset_id_vely, error)
        call h5dcreate_f(file_id, "velz", H5T_NATIVE_DOUBLE, dspace_id, dset_id_velz, error)
#ifdef MHD
        call h5dcreate_f(file_id, "magx", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magx, error)
        call h5dcreate_f(file_id, "magy", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magy, error)
        call h5dcreate_f(file_id, "magz", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magz, error)
#endif !end of ifdef MHD
        call h5dcreate_f(file_id, "pres", H5T_NATIVE_DOUBLE, dspace_id, dset_id_pres, error)
        call h5dcreate_f(file_id, "eint", H5T_NATIVE_DOUBLE, dspace_id, dset_id_eint, error)
        call h5dcreate_f(file_id, "gama", H5T_NATIVE_DOUBLE, dspace_id, dset_id_gama, error)
        call h5dcreate_f(file_id, "ener", H5T_NATIVE_DOUBLE, dspace_id, dset_id_ener, error)
#ifdef MHD
        call h5dcreate_f(file_id, "magp", H5T_NATIVE_DOUBLE, dspace_id, dset_id_magp, error)
        call h5dcreate_f(file_id, "divb", H5T_NATIVE_DOUBLE, dspace_id, dset_id_divb, error)
#endif !end of ifdef MHD

        ! This gets ride of the warning:
        ! "Fortran runtime warning: An array temporary was created"
        ! because we are making the copies of the arrays explicitly
        ! before passing them to the hdf5 subroutine.
        ! However, this does not remove the undesired (and presumabley slow)
        ! behavior that the warning is warning us about.

        dens = grid_V(DENS_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        velx = grid_V(VELX_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        vely = grid_V(VELY_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        velz = grid_V(VELZ_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
#ifdef MHD
        magx = grid_V(MAGX_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        magy = grid_V(MAGY_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        magz = grid_V(MAGZ_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
#endif !end of ifdef MHD
        pres = grid_V(PRES_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        eint = grid_V(EINT_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))
        gama = grid_V(GAMA_VAR, grid_strtIdx(XDIR):grid_stopIdx(XDIR), grid_strtIdx(YDIR):grid_stopIdx(YDIR), grid_strtIdx(ZDIR):grid_stopIdx(ZDIR))

        ! derived quantities
        ener = dens*(velx**2+vely**2+velz**2)/2 + dens*eint ! total energy (hydro)

#ifdef MHD
        ! derived MHD quantities
        magp = (magx**2 + magy**2 + magz**2)/2 ! magnetic pressure
        ener = ener + magp ! total energy (MHD)
        divb = 0.0 ! divergence of B
        do i_dim=1,ndim
            offsets = 0
            offsets(i_dim) = 1
            ii = offsets(XDIR)
            jj = offsets(YDIR)
            kk = offsets(ZDIR)
            divb = divb + 1/(2*grid_dl(i_dim))*&
                    (grid_V(MAGX_VAR+i_dim-1,&
                          grid_strtIdx(XDIR)+ii:grid_stopIdx(XDIR)+ii,&
                          grid_strtIdx(YDIR)+jj:grid_stopIdx(YDIR)+jj,&
                          grid_strtIdx(ZDIR)+kk:grid_stopIdx(ZDIR)+kk)&
                   - grid_V(MAGX_VAR+i_dim-1,&
                          grid_strtIdx(XDIR)-ii:grid_stopIdx(XDIR)-ii,&
                          grid_strtIdx(YDIR)-jj:grid_stopIdx(YDIR)-jj,&
                          grid_strtIdx(ZDIR)-kk:grid_stopIdx(ZDIR)-kk))
        end do
#endif !end of ifdef MHD

        !*****! write data to datasets
        call h5dwrite_f(dset_id_dens, H5T_NATIVE_DOUBLE, dens, data_dims_3d, error)
        call h5dwrite_f(dset_id_velx, H5T_NATIVE_DOUBLE, velx, data_dims_3d, error)
        call h5dwrite_f(dset_id_vely, H5T_NATIVE_DOUBLE, vely, data_dims_3d, error)
        call h5dwrite_f(dset_id_velz, H5T_NATIVE_DOUBLE, velz, data_dims_3d, error)
#ifdef MHD
        call h5dwrite_f(dset_id_magx, H5T_NATIVE_DOUBLE, magx, data_dims_3d, error)
        call h5dwrite_f(dset_id_magy, H5T_NATIVE_DOUBLE, magy, data_dims_3d, error)
        call h5dwrite_f(dset_id_magz, H5T_NATIVE_DOUBLE, magz, data_dims_3d, error)
#endif !end of ifdef MHD
        call h5dwrite_f(dset_id_pres, H5T_NATIVE_DOUBLE, pres, data_dims_3d, error)
        call h5dwrite_f(dset_id_eint, H5T_NATIVE_DOUBLE, eint, data_dims_3d, error)
        call h5dwrite_f(dset_id_gama, H5T_NATIVE_DOUBLE, gama, data_dims_3d, error)
        call h5dwrite_f(dset_id_ener, H5T_NATIVE_DOUBLE, ener, data_dims_3d, error)
#ifdef MHD
        call h5dwrite_f(dset_id_magp, H5T_NATIVE_DOUBLE, magp, data_dims_3d, error)
        call h5dwrite_f(dset_id_divb, H5T_NATIVE_DOUBLE, divb, data_dims_3d, error)
#endif !end of ifdef MHD

        !////! close datasets for primitive variables
        call h5dclose_f(dset_id_dens, error)
        call h5dclose_f(dset_id_velx, error)
        call h5dclose_f(dset_id_vely, error)
        call h5dclose_f(dset_id_velz, error)
#ifdef MHD
        call h5dclose_f(dset_id_magx, error)
        call h5dclose_f(dset_id_magy, error)
        call h5dclose_f(dset_id_magz, error)
#endif !end of ifdef MHD
        call h5dclose_f(dset_id_pres, error)
        call h5dclose_f(dset_id_eint, error)
        call h5dclose_f(dset_id_gama, error)
        call h5dclose_f(dset_id_ener, error)
#ifdef MHD
        call h5dclose_f(dset_id_magp, error)
        call h5dclose_f(dset_id_divb, error)
#endif !end of ifdef MHD

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
        call h5dcreate_f(file_id, "zmin", H5T_NATIVE_DOUBLE, dspace_id, dset_id_zmin, error)
        call h5dcreate_f(file_id, "zmax", H5T_NATIVE_DOUBLE, dspace_id, dset_id_zmax, error)
        call h5dcreate_f(file_id, "nx", H5T_NATIVE_INTEGER, dspace_id, dset_id_nx, error)
        call h5dcreate_f(file_id, "ny", H5T_NATIVE_INTEGER, dspace_id, dset_id_ny, error)
        call h5dcreate_f(file_id, "nz", H5T_NATIVE_INTEGER, dspace_id, dset_id_nz, error)

        !*****! write datasets
        call h5dwrite_f(dset_id_xmin, H5T_NATIVE_DOUBLE, grid_beg(XDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_xmax, H5T_NATIVE_DOUBLE, grid_end(XDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_ymin, H5T_NATIVE_DOUBLE, grid_beg(YDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_ymax, H5T_NATIVE_DOUBLE, grid_end(YDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_zmin, H5T_NATIVE_DOUBLE, grid_beg(ZDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_zmax, H5T_NATIVE_DOUBLE, grid_end(ZDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_nx, H5T_NATIVE_INTEGER, grid_N(XDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_ny, H5T_NATIVE_INTEGER, grid_N(YDIR), data_dims_1d, error)
        call h5dwrite_f(dset_id_nz, H5T_NATIVE_INTEGER, grid_N(ZDIR), data_dims_1d, error)

        !////! close datasets for xmin, xmax, ymin, ymax, nx, ny
        call h5dclose_f(dset_id_xmin, error)
        call h5dclose_f(dset_id_xmax, error)
        call h5dclose_f(dset_id_ymin, error)
        call h5dclose_f(dset_id_ymax, error)
        call h5dclose_f(dset_id_zmin, error)
        call h5dclose_f(dset_id_zmax, error)
        call h5dclose_f(dset_id_nx, error)
        call h5dclose_f(dset_id_ny, error)
        call h5dclose_f(dset_id_nz, error)

        !:::! close dataspace for grid parameters
        call h5sclose_f(dspace_id, error)

        !***! flush all data to disk before closing (helps prevent corruption on interrupt)
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        !++! close the file
        call h5fclose_f(file_id, error)
        
        !=! close hdf5 interface
        call h5close_f(error)
        
    end subroutine output_writeHdf5

end module output
