module diabloH5
    use HDF5
    use MPI
    use grid, only: nx, ny, nz, n_th, rnx, rny, rnz
    use param, only: time, save_flow_int, th, cth
    use decomp_2d, only: xstart, xsize, xend
    use decomp_2d_fft
    implicit none

    private

    public :: write_flow_field
    
contains

!> Subroutine to write the current flow field state out to
!! restart_files/out_xxx.h5
subroutine write_flow_field(final)
    !> Flag determining whether we are writing out flow at end time
    logical, intent(in) :: final

    !> File identifier
    integer(hid_t) :: file_id
    !> Dataset identifier
    integer(hid_t) :: dset_id
    !> Property list identifiers
    integer(hid_t) :: plist_id
    !> Memory space identifiers
    integer(hid_t) :: filespace_id, memspace_id
    !> Dimensions in memory & file of variables
    integer(hsize_t), dimension(3) :: dimsm, dimsf
    !> Starting index offset for local array
    integer(hsize_t), dimension(3) :: offset_f, offset_m
    !> Hyperslab parameters
    integer(hsize_t), dimension(3) :: count, stride
    !> String describing index of output
    character(len= 3) :: frame
    !> Output file name
    character(len=55) :: fname
    !> Dataset name
    character(len=10) :: dname
    !> Number of dimensions
    integer, parameter :: ndims = 3
    !> Error flag for HDF5 functions
    integer :: error
    integer :: n

    dimsm(:) = xsize(:)
    dimsf(:) = [nx, ny, nz]

    offset_f(:) = xstart(:) - 1
    offset_m(:) = 0
    count(:) = 1
    stride(:) = 1


    if (final) then
        fname = 'end.h5'
    else
        write(frame, "(i3.3)") nint(time/save_flow_int)
        fname='restart_files/out'//frame//'.h5'
    end if

    !> Initialize interface
    call h5open_f(error)

    !> Setup file access property list with parallel I/O
    call h5pcreate_f(h5p_file_access_f, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, &
                            mpi_info_null, error)
    !> Create the file
    call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                        file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    !> Convert each variable into physical space
    do n=1,n_th
        call decomp_2d_fft_3d(cth(:,:,:,n), th(:,:,:,n))
    end do
    th = th/rnx/rny/rnz

    ! Missing: create prop list for chunked dataset creation

    !> Create the dataspaces for the variables
    call h5screate_simple_f(ndims, dimsf, filespace_id, error)
    call h5screate_simple_f(ndims, dimsm, memspace_id, error)

    !> Create the property list for the parallel dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)

    do n=1,n_th
        write(frame,'(i1.1)') n
        dname="th"//trim(frame)
        !> Create the dataset item in the file
        call h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, &
                            filespace_id, dset_id, error)
        !> Select hyperslabs in file and memory spaces
        call h5sselect_hyperslab_f(filespace_id, h5s_select_set_f, &
                            offset_f, count, error, stride, dimsm)
        call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                            offset_m, count, error, stride, dimsm)
        !> Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                    th(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3),n), &
                    dimsm, error, file_space_id = filespace_id, mem_space_id = memspace_id)
        !> Close dataset
        call h5dclose_f(dset_id, error)
    end do

    !> Close dataspaces
    call h5sclose_f(filespace_id, error)
    call h5sclose_f(memspace_id, error)

    !> Close property list
    call h5pclose_f(plist_id, error)

    !> Close file and HDF5 interface
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    !> Convert variables back to spectral space
    do n=1,n_th
        call decomp_2d_fft_3d(th(:,:,:,n), cth(:,:,:,n))
    end do

end subroutine write_flow_field
    
end module diabloH5