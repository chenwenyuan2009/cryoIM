subroutine calc_matrix_2f
implicit none
    integer i1,j1
    real r(3,3,60)


    call icos2f_operation_matrix(r,60)
    do i1=1,60
        do j1=1,3
            print*,r(j1,:,i1)
        end do
        write(*,*)
    end do
    return
end


