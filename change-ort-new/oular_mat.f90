subroutine oular_mat(RotMat,theta0,phi0,omega0)
    implicit none
    real theta0,phi0,omega0,deg2rad
    real RotMat(3,3),phi,theta,omega,c1,c2,c3,s1,s2,s3
    parameter(deg2rad=3.1415926/180.0)

    phi  =phi0*deg2rad
    theta=theta0*deg2rad
    omega=omega0*deg2rad
    c1=cos(phi)
    c2=cos(theta)
    c3=cos(omega)
    s1=sin(phi)
    s2=sin(theta)
    s3=sin(omega)
    RotMat(1,1)=c1*c2*c3-s1*s3; RotMat(2,1)=-c3*s1-c1*c2*s3; RotMat(3,1)=c1*s2
    RotMat(1,2)=c1*s3+c2*c3*s1; RotMat(2,2)= c1*c3-c2*s1*s3; RotMat(3,2)=s1*s2
    RotMat(1,3)=-c3*s2;         RotMat(2,3)= s2*s3;          RotMat(3,3)=c2
    return
end subroutine
