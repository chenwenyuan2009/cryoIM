subroutine euler_mat_vertex(theta0,phi0,omega0,RotMat)
    use common_icosproc
    implicit none
    real theta0,phi0,omega0,deg2rad
    real RotMat(3,3),phi1,theta1,omega1,c1,c2,c3,s1,s2,s3
    deg2rad=atan2(1.0,1.0)/45.0
    phi1  =phi0*deg2rad
    theta1=theta0*deg2rad
    omega1=omega0*deg2rad
    c1=cos(phi1)
    c2=cos(theta1)
    c3=cos(omega1)
    s1=sin(phi1)
    s2=sin(theta1)
    s3=sin(omega1)
    RotMat(1,1)=c1*c2*c3-s1*s3; RotMat(2,1)=-c3*s1-c1*c2*s3; RotMat(3,1)=c1*s2
    RotMat(1,2)=c1*s3+c2*c3*s1; RotMat(2,2)= c1*c3-c2*s1*s3; RotMat(3,2)=s1*s2
    RotMat(1,3)=-c3*s2;         RotMat(2,3)= s2*s3;          RotMat(3,3)=c2
    return
end subroutine
