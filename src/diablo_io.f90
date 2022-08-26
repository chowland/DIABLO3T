module IO
    use grid
    use param
    use fft
    use InitialConditions
    use timestepper, only: compute_initial_pressure
    use decomp_2d, only : xstart, xend, nrank
    implicit none
    private

    public :: initialize, finalize
    
contains

subroutine initialize

    call read_input_file
    call create_grid_per
    call init_fft
    call init_vars
    call SetTemperatureIC
    call SetVelocityIC
    call compute_initial_pressure

end subroutine initialize

subroutine finalize

    call destroy_vars
    call finalize_fft

end subroutine finalize

subroutine read_input_file
    real :: current_version = 3.0
    real :: version
    integer :: infile, n
    logical :: exists, reset_time, movie
    character(len=35) :: flavor ! Redundant

    open(newunit=infile, file='input.dat', form='formatted', status='old', action='read')
    read(infile,*)
    read(infile,*)
    read(infile,*)
    read(infile,*)
    read(infile,*) flavor,   version
    if (abs(version - current_version)  > 1e-6) stop 'wrong input data format.'
    read(infile,*)
    read(infile,*) nu, lx, ly, lz
    read(infile,*)
    read(infile,*) create_new_flow
    read(infile,*)
    read(infile,*) n_time_steps, time_limit, delta_t, reset_time, variable_dt, CFL, update_dt
    read(infile,*)
    read(infile,*) verbosity, save_flow_int, save_stats_int, movie
    read(infile,*)
    read(infile,*) !nx_mov, nx_mov_th, ny_mov, ny_mov_th, nz_mov, nz_mov_th
    read(infile,*)
    ! Read in the parameters for the N_TH scalars
    do n = 1,n_th
        read(infile,*)
        read(infile,*) !create_new_th(n)
        read(infile,*)
        read(infile,*) Ri_tau(n), Pr(n)!, reaction(n)
    end do
    close(infile)

    if (nrank==0) write(*,*) 'n_time_steps: ',n_time_steps

    Lx = 2.0_dp*pi*Lx
    Ly = 2.0_dp*pi*Ly
    Lz = 2.0_dp*pi*Lz

    do n=1,n_th
        kappa(n) = nu/Pr(n)
    end do

    dtmax = delta_t

end subroutine read_input_file

end module IO