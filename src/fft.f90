module fft
    use grid
    use param, only : pi, ci, Lx, Ly, Lz
    use decomp_2d
    use decomp_2d_fft
    implicit none
    private

    public :: init_fft, finalize_fft
    real(dp), dimension(1:nx/2+1), public :: kx, kx2
    complex(dp), dimension(1:nx/2+1), public :: cikx
    real(dp), dimension(1:ny), public :: ky, ky2
    complex(dp), dimension(1:ny), public :: ciky
    real(dp), dimension(1:nz), public :: kz, kz2
    complex(dp), dimension(1:nz), public :: cikz
    
contains

subroutine init_fft
    integer :: ierror, i, j, k, nkx, nky, nkz

    call MPI_INIT(ierror)
    call decomp_2d_init(nx,ny,nz,p_row,p_col)
    call decomp_2d_fft_init
    ! Collect start/final indices for spectral space on each process
    call decomp_2d_fft_get_size(cstart,cend,csize)

    nkx = nx/3
    kx(:) = 0.0_dp
    kx2(:) = 0.0_dp
    cikx(:) = 0.0_dp
    do i=0,nkx
        kx(i+1) = i*2.0_dp*pi/Lx
        kx2(i+1) = kx(i+1)**2
        cikx(i+1) = ci*kx(i+1)
    end do
    ! Leave spectral vars equal to zero above cutoff for dealiasing

    nky = ny/3
    ky(:) = 0.0_dp
    ky2(:) = 0.0_dp
    ciky(:) = 0.0_dp
    do j=0,nky
        ky(j+1) = j*2.0_dp*pi/Ly
        ky2(j+1) = ky(j+1)**2
        ciky(j+1) = ci*ky(j+1)
        if (j > 1) then
            ky(ny+1-j) = -ky(j+1)
            ky2(ny+1-j) = ky(j+1)**2
            ciky(ny+1-j) = -ci*ky(j+1)
        end if
    end do

    nkz = nz/3
    kz(:) = 0.0_dp
    kz2(:) = 0.0_dp
    cikz(:) = 0.0_dp
    do k=0,nkz
        kz(k+1) = k*2.0_dp*pi/Lz
        kz2(k+1) = kz(k+1)**2
        cikz(k+1) = ci*kz(k+1)
        if (k > 1) then
            kz(nz+1-k) = -kz(k+1)
            kz2(nz+1-k) = kz(k+1)**2
            cikz(nz+1-k) = -ci*kz(k+1)
        end if
    end do

    if (nrank==0) write(*,*) '2DECOMP FFT package initialized'

end subroutine init_fft

subroutine finalize_fft
    integer :: ierror

    call decomp_2d_fft_finalize
    call decomp_2d_finalize
    call MPI_FINALIZE(ierror)

end subroutine finalize_fft
    
end module fft