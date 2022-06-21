!> Module containing all the global parameters and variables
!> used by the flow solver
module param
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use grid
    use decomp_2d, only : xstart, xend
    implicit none

    real(dp) :: nu, Lx, Ly, Lz, delta_t, CFL
    real(dp) :: save_flow_int, save_stats_int
    real(dp) :: time, time_limit
    real(dp), dimension(1:n_th) :: Ri_tau, Pr, kappa
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), dimension(1:3), parameter :: h_bar = [8.0_dp/15.0_dp, 2.0_dp/15.0_dp, 1.0_dp/3.0_dp]
    real(dp), dimension(1:3), parameter :: beta_bar = [1.0_dp, 25.0_dp/8.0_dp, 9.0_dp/4.0_dp]
    real(dp), dimension(1:3), parameter :: zeta_bar = [0.0_dp, -17.0_dp/8.0_dp, -5.0_dp/4.0_dp]

    complex(dp), parameter :: ci = cmplx(0.0_dp, 1.0_dp, kind=dp)

    integer :: n_time_steps, verbosity, update_dt

    logical :: create_new_flow, variable_dt

    real(dp), allocatable :: th(:,:,:,:), sth(:,:,:)
    complex(dp), allocatable :: cth(:,:,:,:), cfth(:,:,:,:), crth(:,:,:,:), csth(:,:,:)
    
contains

!> Allocate memory for the state variables
!> in both physical and spectral space
!> and set values equal to zero
!> N.B. Run AFTER calling init_fft
subroutine init_vars
    integer :: i, j, k, n

    allocate(  th(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:n_th))
    allocate( sth(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate( cth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate(crth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate(cfth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate(csth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    do n=1,n_th
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    th(i,j,k,n) = 0.0_dp
                    sth(i,j,k) = 0.0_dp
                end do
            end do
        end do
    end do

    do n=1,n_th
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    cth(i,j,k,n) = 0.0_dp
                    csth(i,j,k) = 0.0_dp
                    crth(i,j,k,n) = 0.0_dp
                    cfth(i,j,k,n) = 0.0_dp
                end do
            end do
        end do
    end do

end subroutine init_vars

subroutine create_grid_per

    call create_uniform_grid(gx, Lx, nx)
    call create_uniform_grid(gy, Ly, ny)
    call create_uniform_grid(gz, Lz, nz)

end subroutine create_grid_per

subroutine destroy_vars
    deallocate(th,sth)
    deallocate(cth)
    deallocate(crth)
    deallocate(cfth)
end subroutine destroy_vars

end module param