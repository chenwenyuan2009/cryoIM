subroutine icos2f_sym_matrix(r,N)
    implicit none
    integer N,i,j
    real u,v,w,deg2rad,a5
    real r(3,3,N),rz5(3,3),ryv(3,3),rz15(3,3),rx2(3,3),ry2(3,3),r111_3(3,3)

    deg2rad=4.0*atan2(1.0,1.0)/180.0
    a5=72.0*deg2rad
    u=atan(1.0/tan(60.0*deg2rad)/tan(36.0*deg2rad)) !angle between 3f and 5f
    v=acos(cos(60.0*deg2rad)/sin(36.0*deg2rad))     !angle between 2f and 5f
    w=acos(cos(36.0*deg2rad)/sin(60.0*deg2rad))     !angle between 2f and 3f

    rz5(1,1)=cos(a5);rz5(1,2)=-sin(a5);rz5(1,3)= .0
    rz5(2,1)=sin(a5);rz5(2,2)= cos(a5);rz5(2,3)= .0
    rz5(3,1)=.0     ;rz5(3,2)=.0      ;rz5(3,3)=1.0

    ryv(1,1)= cos(v);ryv(1,2)=.0;ryv(1,3)=sin(v)
    ryv(2,1)=.0     ;ryv(2,2)=1.;ryv(2,3)=.0
    ryv(3,1)=-sin(v);ryv(3,2)=.0;ryv(3,3)=cos(v)

    rz5=matmul(ryv,matmul(rz5,transpose(ryv)))

    rx2(1,1)=1.0;rx2(1,2)=  .0;rx2(1,3)=  .0
    rx2(2,1)= .0;rx2(2,2)=-1.0;rx2(2,3)=  .0
    rx2(3,1)= .0;rx2(3,2)=  .0;rx2(3,3)=-1.0

    ry2(1,1)=-1.0;ry2(1,2)=  .0;ry2(1,3)=  .0
    ry2(2,1)=  .0;ry2(2,2)= 1.0;ry2(2,3)=  .0
    ry2(3,1)=  .0;ry2(3,2)=  .0;ry2(3,3)=-1.0

    r111_3(1,1)= .0;r111_3(1,2)= .0;r111_3(1,3)=1.0
    r111_3(2,1)=1.0;r111_3(2,2)= .0;r111_3(2,3)= .0
    r111_3(3,1)= .0;r111_3(3,2)=1.0;r111_3(3,3)= .0

    r=.0
    r(1,1,1)=1.0;r(2,2,1)=1.0;r(3,3,1)=1.0
    do i=2,5
        r(:,:,i)=matmul(r(:,:,i-1),rz5)
    end do

    do i=6,10
        r(:,:,i)=matmul(rx2,r(:,:,i-5))
    end do

    do i=11,20
        r(:,:,i)=matmul(ry2,r(:,:,i-10))
    end do
    
    do i=21,40
        r(:,:,i)=matmul(r111_3,r(:,:,i-20))
    end do

    do i=41,60
        r(:,:,i)=matmul(r111_3,r(:,:,i-20))
    end do

    return
end
        
        