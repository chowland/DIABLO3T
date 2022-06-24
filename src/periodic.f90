module timestepper
    use param
    use fft
    use decomp_2d
    use decomp_2d_fft
    implicit none
    private

    public :: rk_per_1, compute_initial_pressure
    
contains

subroutine rk_per_1(rk_step)
    integer, intent(in) :: rk_step
    
    integer :: i, j, k, n
    real(dp) :: temp1, temp2, temp3, temp4, temp5

    temp4 = h_bar(rk_step)*delta_t
    temp1 = nu*temp4*0.5_dp
    temp2 = beta_bar(rk_step)*temp4
    temp3 = zeta_bar(rk_step)*temp4

    ! Explicit part of viscous term
    do n=1,3
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp5 = 1.0_dp - temp1*(kx2(i) + ky2(j) + kz2(k))
                    crhs(i,j,k,n) = temp5*cu(i,j,k,n)
                end do
            end do
        end do
    end do
    ! Add pressure gradient (explicit Euler implementation for each substep)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                crhs(i,j,k,1) = crhs(i,j,k,1) - temp4*cikx(i)*cp(i,j,k)
                crhs(i,j,k,2) = crhs(i,j,k,2) - temp4*ciky(j)*cp(i,j,k)
                crhs(i,j,k,3) = crhs(i,j,k,3) - temp4*cikz(k)*cp(i,j,k)
            end do
        end do
    end do

    ! Explicit part of diffusive term for scalars
    do n=1,n_th
        temp1 = 0.5_dp*kappa(n)*temp4
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp5 = 1.0_dp - temp1*(kx2(i) + ky2(j) + kz2(k))
                    crth(i,j,k,n) = temp5*cth(i,j,k,n)
                end do
            end do
        end do
    end do

    if (rk_step > 1) then
        ! Loop over velocity components
        do n=1,3
            do k=cstart(3),cend(3)
                do j=cstart(2),cend(2)
                    do i=cstart(1),cend(1)
                        crhs(i,j,k,n) = crhs(i,j,k,n) + temp3*cfu(i,j,k,n)
                    end do
                end do
            end do
        end do
        ! Loop over scalar fields
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
    
    ! Advection in momentum equation
    !!! NB DIABLO_PER USES JUST 6 FFTS HERE, TAKING ADVANTAGE OF SYMMETRY
    !!! THIS WILL PROBABLY WANT IMPLEMENTING
    do n=1,3
        ! Reset array to zero (add forcing and buoyancy terms here if using)
        cfu(:,:,:,n) = 0.0_dp
        call advect_field(u(:,:,:,n), cfu(:,:,:,n))
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    crhs(i,j,k,n) = crhs(i,j,k,n) + temp2*cfu(i,j,k,n)
                end do
            end do
        end do
    end do

    ! Advection in scalar equations
    do n=1,n_th
        ! Reset array to zero (add forcing terms here if using)
        cfth(:,:,:,n) = 0.0_dp
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
    temp1 = 0.5_dp*nu*temp4
    do n=1,3
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp5 = 1.0_dp + temp1*(kx2(i) + ky2(j) + kz2(k))
                    cu(i,j,k,n) = crhs(i,j,k,n)/temp5
                end do
            end do
        end do
    end do

    ! Solve the implicit system for the scalar at next timestep
    do n=1,n_th
        temp1 = 0.5_dp*kappa(n)*temp4
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    temp5 = 1 + temp1*(kx2(i) + ky2(j) + kz2(k))
                    cth(i,j,k,n) = crth(i,j,k,n)/temp5
                end do
            end do
        end do
    end do

    ! Compute second fractional step...
    ! Make cu divergence free, and save pressure step in cs1
    call remove_divergence
    ! ! ! Update pressure
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cp(i,j,k) = cp(i,j,k) + cs1(i,j,k)/temp4
            end do
        end do
    end do

    ! Update the velocity fields in physical space
    do n=1,n_th
        call decomp_2d_fft_3d(cu(:,:,:,n), u(:,:,:,n))
    end do
    u = u/rnx/rny/rnz

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

    do n=1,3
        ! x-advection
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    s1(i,j,k) = u(i,j,k,n)*var(i,j,k)
                end do
            end do
        end do
        call decomp_2d_fft_3d(s1, cs1)
        ! Dealias the highest modes
        do k=cstart(3),cend(3)
            do j=cstart(2),cend(2)
                do i=cstart(1),cend(1)
                    if (alias(i,j,k)) cs1(i,j,k) = 0.0_dp
                end do
            end do
        end do
        if (n==1) then
            do k=cstart(3),cend(3)
                do j=cstart(2),cend(2)
                    do i=cstart(1),cend(1)
                        cfvar(i,j,k) = cfvar(i,j,k) - cikx(i)*cs1(i,j,k)
                    end do
                end do
            end do
        else if (n==2) then
            do k=cstart(3),cend(3)
                do j=cstart(2),cend(2)
                    do i=cstart(1),cend(1)
                        cfvar(i,j,k) = cfvar(i,j,k) - ciky(j)*cs1(i,j,k)
                    end do
                end do
            end do
        else if (n==3) then
            do k=cstart(3),cend(3)
                do j=cstart(2),cend(2)
                    do i=cstart(1),cend(1)
                        cfvar(i,j,k) = cfvar(i,j,k) - cikz(k)*cs1(i,j,k)
                    end do
                end do
            end do
        end if
    end do
    
end subroutine advect_field
    
subroutine remove_divergence
    integer :: i, j, k
    real(dp) :: temp

    cs1(:,:,:) = 0.0_dp

    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                temp = -(kx2(i) + ky2(j) + kz2(k) + 1e-14)
                cs1(i,j,k) = (cikx(i)*cu(i,j,k,1) + ciky(j)*cu(i,j,k,2) &
                                + cikz(k)*cu(i,j,k,3))/temp
                cu(i,j,k,1) = cu(i,j,k,1) - cikx(i)*cs1(i,j,k)
                cu(i,j,k,2) = cu(i,j,k,2) - ciky(j)*cs1(i,j,k)
                cu(i,j,k,3) = cu(i,j,k,3) - cikz(k)*cs1(i,j,k)
            end do
        end do
    end do

end subroutine remove_divergence

subroutine compute_initial_pressure
    integer :: i, j, k, n
    real(dp), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),2) :: s2

    ! Use s1 as working variable for pressure in physical space
    s1(:,:,:) = 0.0_dp

    ! (du/dx) & (dv/dy) [& implicitly (dw/dz)]
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfu(i,j,k,1) = cikx(i)*cu(i,j,k,1)
                cfu(i,j,k,2) = ciky(j)*cu(i,j,k,2)
            end do
        end do
    end do
    do n=1,2
        call decomp_2d_fft_3d(cfu(:,:,:,n), s2(:,:,:,n))
    end do
    s2 = s2/rnx/rny/rnz
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = 2*(s2(i,j,k,1)*s2(i,j,k,1) &
                        + s2(i,j,k,1)*s2(i,j,k,2) + s2(i,j,k,2)*s2(i,j,k,2))
            end do
        end do
    end do
    ! (du/dy) & (dv/dx)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfu(i,j,k,1) = ciky(j)*cu(i,j,k,1)
                cfu(i,j,k,2) = cikx(i)*cu(i,j,k,2)
            end do
        end do
    end do
    do n=1,2
        call decomp_2d_fft_3d(cfu(:,:,:,n), s2(:,:,:,n))
    end do
    s2 = s2/rnx/rny/rnz
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = s1(i,j,k) + 2*s2(i,j,k,1)*s2(i,j,k,2)
            end do
        end do
    end do
    ! (du/dz) & (dw/dx)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfu(i,j,k,1) = cikz(k)*cu(i,j,k,1)
                cfu(i,j,k,2) = cikx(i)*cu(i,j,k,3)
            end do
        end do
    end do
    do n=1,2
        call decomp_2d_fft_3d(cfu(:,:,:,n), s2(:,:,:,n))
    end do
    s2 = s2/rnx/rny/rnz
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = s1(i,j,k) + 2*s2(i,j,k,1)*s2(i,j,k,2)
            end do
        end do
    end do
    ! (dv/dz) & (dw/dy)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                cfu(i,j,k,1) = cikz(k)*cu(i,j,k,2)
                cfu(i,j,k,2) = ciky(j)*cu(i,j,k,3)
            end do
        end do
    end do
    do n=1,2
        call decomp_2d_fft_3d(cfu(:,:,:,n), s2(:,:,:,n))
    end do
    s2 = s2/rnx/rny/rnz
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = s1(i,j,k) + 2*s2(i,j,k,1)*s2(i,j,k,2)
            end do
        end do
    end do

    call decomp_2d_fft_3d(s1, cp)
    do k=cstart(3),cend(3)
        do j=cstart(2),cend(2)
            do i=cstart(1),cend(1)
                if (alias(i,j,k)) cp(i,j,k) = 0.0
                cp(i,j,k) = cp(i,j,k)/(kx2(i) + ky2(j) + kz2(k) + 1e-14)
            end do
        end do
    end do

end subroutine compute_initial_pressure

end module timestepper