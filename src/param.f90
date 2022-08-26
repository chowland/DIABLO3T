!> Module containing all the global parameters and variables
!> used by the flow solver
module param
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use grid
    use decomp_2d, only : xstart, xend
    implicit none

    real(dp) :: nu, Lx, Ly, Lz, delta_t, CFL, dtmax
    real(dp) :: save_flow_int, save_stats_int
    real(dp) :: time, time_limit
    real(dp), dimension(1:n_th) :: Ri_tau, Pr, kappa
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) !< This is pi
    real(dp), dimension(1:3), parameter :: h_bar = [8.0_dp/15.0_dp, 2.0_dp/15.0_dp, 1.0_dp/3.0_dp]
    real(dp), dimension(1:3), parameter :: beta_bar = [1.0_dp, 25.0_dp/8.0_dp, 9.0_dp/4.0_dp]
    real(dp), dimension(1:3), parameter :: zeta_bar = [0.0_dp, -17.0_dp/8.0_dp, -5.0_dp/4.0_dp]

    complex(dp), parameter :: ci = cmplx(0.0_dp, 1.0_dp, kind=dp)
        !< The imaginary unit

    integer :: n_time_steps, verbosity, update_dt

    logical :: create_new_flow, variable_dt

    real(dp), allocatable :: u(:,:,:,:)  !< velocity field
    real(dp), allocatable :: th(:,:,:,:)    !< Scalar field values
    complex(dp), allocatable :: cth(:,:,:,:)    !< Scalar fields in spectral space
    complex(dp), allocatable :: cu(:,:,:,:) !< velocity field in spectral space
    complex(dp), allocatable :: cp(:,:,:)   !< pressure field in spectral space

    real(dp), allocatable :: s1(:,:,:)
        !< Temporary variable in physical space used in timestepper
        !! for nonlinear term calculation
    complex(dp), allocatable :: cs1(:,:,:)
        !< Temporary spectral space array
    complex(dp), allocatable :: cfu(:,:,:,:)
        !< Array storing nonlinear terms for RK timestepping
        !! Stores data from previous sub-step for use
        !! so DO NOT OVERWRITE outside of timestepper
    complex(dp), allocatable :: cfth(:,:,:,:)
        !< Array storing nonlinear terms for RK timestepping
        !! Stores data from previous sub-step for use
        !! so DO NOT OVERWRITE outside of timestepper
    complex(dp), allocatable :: crhs(:,:,:,:)
    complex(dp), allocatable :: crth(:,:,:,:)
        !< Array used for rhs of scalar evolution equation
    
contains

!> Allocate memory for the state variables
!> in both physical and spectral space
!> and set values equal to zero
!> N.B. Run AFTER calling init_fft
subroutine init_vars
    ! integer :: i, j, k, n

    !> Allocate memory for velocity arrays
    allocate(   u(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:3))
    allocate(  cu(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:3))
    allocate(  cp(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    allocate(crhs(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:3))
    allocate( cfu(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:3))

    allocate(  th(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:n_th))
    allocate(  s1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate( cth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate(crth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate(cfth(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),1:n_th))
    allocate( cs1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

end subroutine init_vars

subroutine create_grid_per

    call create_uniform_grid(gx, Lx, nx)
    call create_uniform_grid(gy, Ly, ny)
    call create_uniform_grid(gz, Lz, nz)

end subroutine create_grid_per

subroutine destroy_vars
    deallocate(u)
    deallocate(cu, cp)
    deallocate(s1, cs1)
    deallocate(th, cth)
    deallocate(crhs, cfu)
    deallocate(crth, cfth)
end subroutine destroy_vars

end module param