subroutine icos_matrix_I3(r,N)
implicit none
    integer N,i1,j1
    real r(3,3,N),v,u,w,deg2rad,rot0(3,3),rad2deg

    deg2rad=4.0*atan2(1.0,1.0)/180.0
    rad2deg=180.0/4.0/atan2(1.0,1.0)
    v=acos(cos(60.0*deg2rad)/sin(36.0*deg2rad))     !angle between 2f and 5f
    u=acos(1.0/tan(60.0*deg2rad)/tan(36.0*deg2rad)) !angle between 3f and 5f
    w=acos(cos(36.0*deg2rad)/sin(60.0*deg2rad))     !angle between 2f and 3f
    rot0(1,1)=1.; rot0(1,2)=.0     ; rot0(1,3)=.0
    rot0(2,1)=.0; rot0(2,2)= cos(-w); rot0(2,3)=sin(-w)
    rot0(3,1)=.0; rot0(3,2)=-sin(-w); rot0(3,3)=cos(-w)

    call icos_matrix_I2(r,N)
    do i1=1,N
        r(:,:,i1)=matmul(rot0,matmul(r(:,:,i1),transpose(rot0)))
    end do
    return
end


