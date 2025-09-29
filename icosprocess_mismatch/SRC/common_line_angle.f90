
subroutine  crosslin_angle(theta1,phi1,omega1,theta2,phi2,omega2,alpha1,alpha2)
    implicit none
    integer i
    real rad2deg,deg2rad,theta1,theta2,phi1,phi2,omega1,omega2
    real rot(3,3,60),euler_mat1(3,3),euler_mat2(3,3),alpha1,alpha2
    real z0(3),z1(3),z2(3),z2a(3),vec1(3),vv,vec10(3),vec20(3)

    rad2deg=45.0/atan2(1.0,1.0)
    deg2rad=1.0/rad2deg
    call euler_mat(theta1,phi1,omega1,euler_mat1)
    call euler_mat(theta2,phi2,omega2,euler_mat2)
    z0(1)=.0;z0(2)=.0;z0(3)=1.
    z1=matmul(euler_mat1,z0)  !image
    z2=matmul(euler_mat2,z0)  !template

        vec1(1)=z1(2)*z2(3)-z1(3)*z2(2)
        vec1(2)=z1(3)*z2(1)-z1(1)*z2(3)
        vec1(3)=z1(1)*z2(2)-z1(2)*z2(1)
        vv=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
        if(vv.lt..5) then
            alpha1=1000.0
            alpha2=1000.0
            goto 200
        end if
        vec1 =vec1/vv
        vec10=matmul(transpose(euler_mat1),vec1)
        vec20=matmul(transpose(euler_mat2),vec1)
        alpha1=atan2(vec10(2),vec10(1))*rad2deg
        alpha2=atan2(vec20(2),vec20(1))*rad2deg


200  end subroutine
