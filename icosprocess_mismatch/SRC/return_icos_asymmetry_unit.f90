subroutine return_icos_asymmetry_unit(theta,phi,omega)
    implicit none
    integer i
    real theta,phi,omega,rot(3,3,60),mat(3,3),matoular(3,3)
    real A1,B1,C1,A2,B2,C2,sin31,cos31,sin69,cos69,deg2rad
    real x,y,z,z1,z2

    deg2rad=atan2(1.0,1.0)/45.0
    cos31=cos(31.717*deg2rad);sin31=sin(31.717*deg2rad)
    cos69=cos(69.090*deg2rad);sin69=sin(69.090*deg2rad)
    A1=-sin31*cos69;B1=-cos31*cos69;C1=sin31*sin69
    A2=-sin31*cos69;B2= cos31*cos69;C2=sin31*sin69

    call icos2f_sym_mat(rot)
    call euler_mat(theta,phi,omega,matoular)

    do i=1,60
        mat=matmul(rot(:,:,i),matoular)
        theta=acos(mat(3,3))
        phi  =atan2(mat(2,3),mat(1,3))
        omega=atan2(mat(3,2),-mat(3,1))
        x=sin(theta)*cos(phi)
        y=sin(theta)*sin(phi)
        z=cos(theta)
        z1=-(A1*x+B1*y)/C1
        z2=-(A2*x+B2*y)/C2
        if((z.ge..0).and.(z.le.z1).and.(z.le.z2)) goto 100
    end do
100 theta=theta/deg2rad
    phi  =phi/deg2rad
    omega=omega/deg2rad
    return
end



