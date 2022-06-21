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
    temp1 = nu*temp4*0.5_dp
    temp2 = beta_bar(rk_step)*temp4
    temp3 = zeta_bar(rk_step)*temp4

    ! Explicit part of diffusive term
    do n=1,n_th
        temp1 = 0.5_dp*kappa(n)*temp4
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
    

    !! ADVECT THE TRACER - NONLINEAR TERMS
    do n=1,n_th
        ! Construct nonlinear terms in cfth
        call advect_field(th(:,:,:,n), cfth(:,:,:,n))
        ! Add nonlinear terms to rhs
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    crth(i,j,k,n) = crth(i,j,k,n) + temp2*cfth(i,j,k,n)
                end do
            end do
        end do
    end do

    ! Solve the implicit system for the intermediate field
    do n=1,n_th
        temp1 = 0.5_dp*kappa(n)*temp4
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

    ! Update the scalar fields in physical space at the end of the time step
    do n=1,n_th
        call decomp_2d_fft_3d(cth(:,:,:,n), th(:,:,:,n))
    end do
    th = th/rnx/rny/rnz

end subroutine rk_per_1

!> Subroutine for updating cfth with the nonlinear (advection) terms
!> RHS array crth is updated with the input factor temp * nonlinear terms
subroutine advect_field(var, cfvar)
    !> Field (in physical space) to be advected
    real(dp), dimension(xstart(1):,xstart(2):,xstart(3):), intent(in) :: var
    !> Output field: (minus) divergence of (velocity * var) in spectral space
    complex(dp), dimension(cstart(1):,cstart(2):,cstart(3):), intent(out) :: cfvar

    integer :: i, j, k, n

    cfvar(:,:,:) = 0.0_dp
    ! x-advection
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = u1(i,j,k)*var(i,j,k)
            end do
        end do
    end do
    call decomp_2d_fft_3d(s1, cs1)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfvar(i,j,k) = cfvar(i,j,k) - cikx(i)*cs1(i,j,k)
            end do
        end do
    end do
    ! y-advection
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = u2(i,j,k)*var(i,j,k)
            end do
        end do
    end do
    call decomp_2d_fft_3d(s1, cs1)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfvar(i,j,k) = cfvar(i,j,k) - ciky(j)*cs1(i,j,k)
            end do
        end do
    end do
    ! z-advection
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = u3(i,j,k)*var(i,j,k)
            end do
        end do
    end do
    call decomp_2d_fft_3d(s1, cs1)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfvar(i,j,k) = cfvar(i,j,k) - cikz(k)*cs1(i,j,k)
            end do
        end do
    end do

end subroutine advect_field
    
end module timestepper