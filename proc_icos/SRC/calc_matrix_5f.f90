subroutine calc_matrix_5f
implicit none
    integer i1,j1
    real r(3,3,60),v,deg2rad,rot0(3,3),rad2deg

    deg2rad=4.0*atan2(1.0,1.0)/180.0
    rad2deg=180.0/4.0/atan2(1.0,1.0)
    v=acos(cos(60.0*deg2rad)/sin(36.0*deg2rad))     !angle between 2f and 5f
    rot0(1,1)=cos(v) ;rot0(1,2)=.0;rot0(1,3)=sin(v)
    rot0(2,1)=.0      ;rot0(2,2)=1.;rot0(2,3)=.0
    rot0(3,1)=-sin(v);rot0(3,2)=.0;rot0(3,3)=cos(v)

    call icos2f_operation_matrix(r,60)
    do i1=1,60
        r(:,:,i1)=matmul(rot0,matmul(r(:,:,i1),transpose(rot0)))
        do j1=1,3
            print*,r(j1,:,i1)
        end do
        write(*,*)
    end do
    return
end


