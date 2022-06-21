module grid
    use, intrinsic :: iso_fortran_env, only : dp => real64
    ! use param
    implicit none
    
    !!! USER INPUT !!!
    ! (Hard-coded to improve optimizations)
    ! (Re-compile after any changes)

    ! Grid dimensions
    integer, parameter :: nx = 128  !< Grid size in x
    integer, parameter :: ny = 128  !< Grid size in y
    integer, parameter :: nz = 1    !< Grid size in z
    integer, parameter :: n_th = 2  !< Number of scalar variables
    integer, parameter :: p_row=4, p_col=1  !< Pencil decomposition factors

    !!! END OF USER INPUT !!!

    real(dp) :: gx(1:nx), gy(1:ny), gz(1:nz) !< Grid coordinates
    !< Grid coordinates
    real(dp) :: dx, dy, dz
    ! Real value versions of grid sizes (for FFT normalisation)
    real(dp), parameter :: rnx = real(nx, kind=dp)
    real(dp), parameter :: rny = real(ny, kind=dp)
    real(dp), parameter :: rnz = real(nz, kind=dp)
    ! Indices for each process in spectral space
    ! (initialized in init_fft)
    integer, dimension(3) :: cstart, cend, csize

contains

subroutine create_uniform_grid(gc, L, n)
    !> Grid coordinate vector
    real(dp), intent(out) :: gc(:)
    !> Physical length of grid
    real(dp), intent(in) :: L
    !> Number of grid points
    integer, intent(in) :: n

    integer :: i

    do i=1,n
        gc(i) = (i*L)/n
        ! if (verbosity > 3) write(*,*) 'gx(',i,') = ',gx(i)
    end do
    dx = L/n

end subroutine create_uniform_grid

end module grid