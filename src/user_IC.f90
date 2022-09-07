module InitialConditions
    use param
    use grid
    use fft, only: alias
    use timestepper, only: remove_divergence
    use decomp_2d, only: xstart, xend, nrank
    use decomp_2d_fft

    private
    public :: SetTemperatureIC, SetVelocityIC, SetPressure
contains

subroutine SetTemperatureIC
    integer :: i, j, k, n

    call random_seed()

    if (n_th > 0) then
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    call random_number(rvar)
                    rvar = 2*rvar - 1.0_dp
                    th(i,j,k,1) = 1e-1*rvar
                    ! th(i,j,k,1) = 0.1 * sin(2*pi*(3*gx(i)/Lx + gy(j)/Ly))
                end do
            end do
        end do
    end if

    if (n_th > 1) then
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    call random_number(rvar)
                    rvar = 2*rvar - 1.0_dp
                    th(i,j,k,2) = 1e-1*rvar
                    ! th(i,j,k,2) = 0.1*cos(2*pi*(gx(i)/Lx + gy(j)/Ly))
                end do
            end do
        end do
    end if

    do n=1,n_th
        call decomp_2d_fft_3d(th(:,:,:,n), cth(:,:,:,n))
    end do
end subroutine SetTemperatureIC

subroutine SetVelocityIC
    integer :: i, j, k
    real(dp) :: amp, rvar

    call random_seed()

    amp = 1e-1
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                call random_number(rvar)
                rvar = 2*rvar - 1.0_dp
                ! u(i,j,k,1) = cos(gx(i)) * sin(gy(j)) + amp*rvar
                ! u(i,j,k,1) = sin(gy(j)) + amp*rvar
                u(i,j,k,1) = amp*rvar
                call random_number(rvar)
                rvar = 2*rvar - 1.0_dp
                ! u(i,j,k,2) =-sin(gx(i)) * cos(gy(j)) + amp*rvar
                u(i,j,k,2) = amp*rvar
            end do
        end do
    end do

    u(:,:,:,3) = 0.0_dp
    u(:,:,:,:) = 0.0_dp

    do n=1,3
        call decomp_2d_fft_3d(u(:,:,:,n), cu(:,:,:,n))
    end do

    call remove_divergence

    ! do k=cstart(3),cend(3)
    !     do j=cstart(2),cend(2)
    !         do i=cstart(1),cend(1)
    !             if (alias(i,j,k)) cu(i,j,k,:) = 0.0
    !             if ((i==1) .and. (j==1) .and. (k==1)) cu(i,j,k,:) = 0.0
    !         end do
    !     end do
    ! end do

    ! if (nrank==0) write(*,*) 'cu(1,:,1,1): ', cu(1,:,1,1)
    ! if (nrank==0) write(*,*) 'cu(2,:,1,1): ', cu(2,:,1,1)

    ! if (nrank==0) write(*,*) 'cu(1,:,1,2): ', cu(1,:,1,2)
    ! if (nrank==0) write(*,*) 'cu(2,:,1,2): ', cu(2,:,1,2)

end subroutine SetVelocityIC

end module InitialConditions