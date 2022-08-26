program DIABLO
    use param
    use IO
    use timestepper
    use mpi
    use decomp_2d, only : nrank, nproc
    use diabloH5

    implicit none
    
    integer :: time_step, first_step
    integer :: rk_step, ierror
    logical :: end_flag
    real :: start_time, end_time
    
    end_flag = .false.

    call initialize

    call write_flow_field(.false.)

    start_time = mpi_wtime()

    if (nrank == 0) then
        write(*,*)
        write(*,*) '             ****** Welcome to DIABLO ******'
        write(*,*)
        write(*,*) 'MPI initialized with ', nproc, ' processes'
    end if

    first_step = time_step

    do time_step = time_step+1, time_step+n_time_steps
        if (nrank == 0) write(*,*) 'Now beginning time step = ',time_step,'  dt = ',delta_t

        call update_dt_CFL

        do rk_step = 1,3
            call rk_per_1(rk_step)
        end do

        time = time + delta_t

        if (mod(time, save_flow_int) < delta_t) call write_flow_field(.false.)

        ! Check if exceeding maximum wall-time
        if (nrank == 0) then
            end_time = mpi_wtime()
            if (end_time - start_time > time_limit) then
                write(*,*) 'STOP because of wall-time hit!'
                end_flag = .true.
            end if
        end if
        call mpi_bcast(end_flag, 1, mpi_logical, 0, mpi_comm_world, ierror)
        if (end_flag) exit
    end do

    n_time_steps = time_step - first_step

    if (nrank == 0) then
        write(*,*) 'Elapsed time (sec): ', end_time - start_time
        write(*,*) 'Seconds per iteration: ', (end_time - start_time)/n_time_steps
    end if

    call finalize

end program DIABLO