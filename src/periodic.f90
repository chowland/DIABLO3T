module timestepper
    use param
    use fft
    use decomp_2d
    use decomp_2d_fft
    implicit none
    private

    public :: rk_per_1
    
contains

subroutine rk_per_1(rk_step)
    integer, intent(in) :: rk_step
    
    integer :: i, j, k, n
    real(dp) :: temp1, temp2, temp3, temp4, temp6

    temp4 = h_bar(rk_step)*delta_t
    temp1 = kappa(1)*temp4/2.0_dp
    temp2 = beta_bar(rk_step)*temp4
    temp3 = zeta_bar(rk_step)*temp4

    ! Explicit part of diffusive term
    do n=1,n_th
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp6 = 1.0_dp - temp1*(kx2(i) + ky2(j) + kz2(k))
                    crth(i,j,k,n) = temp6*cth(i,j,k,n)
                end do
            end do
        end do
    end do

    if (rk_step > 1) then
        do n=1,n_th
            do k=cstart(3),cend(3)
                do j=cstart(2),cend(2)
                    do i=cstart(1),cend(1)
                        crth(i,j,k,n) = crth(i,j,k,n) + temp3*cfth(i,j,k,n)
                    end do
                end do
            end do
        end do
    end if

    ! Transform the scalar concentration to physical space
    do n=1,n_th
        call decomp_2d_fft_3d(cth(:,:,:,n), th(:,:,:,n))
    end do
    th = th/rnx/rny/rnz

    !! ADVECT THE TRACER - NONLINEAR TERMS
    do n=1,n_th
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    sth(i,j,k) = sin(2*pi/Ly*gy(j))*th(i,j,k,n)
                end do
            end do
        end do
        call decomp_2d_fft_3d(sth, csth)
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    cfth(i,j,k,n) = cfth(i,j,k,n) - cikx(i)*csth(i,j,k)
                    crth(i,j,k,n) = crth(i,j,k,n) + temp2*cfth(i,j,k,n)
                end do
            end do
        end do
        
    end do

    ! Solve the implicit system for the intermediate field
    do n=1,n_th
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp6 = 1 + temp1*(kx2(i) + ky2(j) + kz2(k))
                    cth(i,j,k,n) = crth(i,j,k,n)/temp6
                end do
            end do
        end do
    end do

    ! Compute second fractional step...


end subroutine rk_per_1
    
end module timestepper