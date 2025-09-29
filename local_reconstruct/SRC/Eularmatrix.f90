subroutine Eularmatrix(theta0,phi0,omega0,r)
    implicit none
    real theta0,theta,phi0,phi,omega0,omega,deg2rad,r(3,3)
    
    deg2rad=4.0*atan2(1.0,1.0)/180.0
    theta=theta0*deg2rad
    phi  =phi0  *deg2rad
    omega=omega0*deg2rad

    r=.0
    r(1,1)=cos(phi)*cos(theta)*cos(omega)-sin(phi)*sin(omega);r(1,2)=-cos(phi)*cos(theta)*sin(omega)-sin(phi)*cos(omega);r(1,3)=cos(phi)*sin(theta)
    r(2,1)=sin(phi)*cos(theta)*cos(omega)+cos(phi)*sin(omega);r(2,2)=-sin(phi)*cos(theta)*sin(omega)+cos(phi)*cos(omega);r(2,3)=sin(phi)*sin(theta)
    r(3,1)=-sin(theta)*cos(omega)                            ;r(3,2)= sin(theta)*sin(omega)                             ;r(3,3)=cos(theta)
    return
end


