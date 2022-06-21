module InitialConditions
    use param
    use grid
    use decomp_2d, only: xstart, xend
    use decomp_2d_fft

    private
    public :: SetTemperatureIC
contains

subroutine SetTemperatureIC
    integer :: i, j, k, n

    do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do i=xstart(1),xend(1)
                th(i,j,k,1) = sin(2*pi*(3*gx(i)/Lx + 4*gy(j)/Ly))
            end do
        end do
    end do

    if (n_th > 1) then
        do k=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do i=xstart(1),xend(1)
                    th(i,j,k,2) = cos(2*pi*(5*gx(i)/Lx + 2*gy(j)/Ly))
                end do
            end do
        end do
    end if

    do n=1,n_th
        call decomp_2d_fft_3d(th(:,:,:,n), cth(:,:,:,n))
    end do
end subroutine SetTemperatureIC

end module InitialConditions