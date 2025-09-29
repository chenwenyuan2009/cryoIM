

subroutine euler_mat(theta0,phi0,omega0,RotMat)
    implicit none
    real RotMat(3,3),theta0,theta,phi0,phi,omega0,omega,c1,c2,c3,s1,s2,s3,deg2rad
    deg2rad=atan2(1.0,1.0)/45.0
    theta=theta0*deg2rad
    phi  =phi0  *deg2rad
    omega=omega0*deg2rad
    c1=cos(phi)
    c2=cos(theta)
    c3=cos(omega)
    s1=sin(phi)
    s2=sin(theta)
    s3=sin(omega)
    RotMat(1,1)=c1*c2*c3-s1*s3; RotMat(1,2)=-c1*c2*s3-c3*s1; RotMat(1,3)=c1*s2
    RotMat(2,1)=c2*c3*s1+c1*s3; RotMat(2,2)=-c2*s1*s3+c1*c3; RotMat(2,3)=s1*s2
    RotMat(3,1)=-c3*s2;         RotMat(3,2)= s2*s3;          RotMat(3,3)=c2

    return
end subroutine

