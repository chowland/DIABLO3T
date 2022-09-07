module forcing
    use grid
    use param, only: s1, cfu
    use decomp_2d, only: xstart, xend
    use decomp_2d_fft
    implicit none
    
contains

!> Subroutine for adding custom forcing to the momentum equation
!! Use the temporary array s1 to define the forcing in physical space
!! and then transform to the desired component of cfu
subroutine add_momentum_forcing
    integer :: i, j, k

    ! sin(2y) forcing for Kolmogorov flow
    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                s1(i,j,k) = 0.16*cos(gy(j))
            end do
        end do
    end do

    call decomp_2d_fft_3d(s1, cfu(:,:,:,1))

end subroutine
    
end module forcing