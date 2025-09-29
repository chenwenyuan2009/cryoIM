subroutine get_randomort
    use common_icosproc
    implicit none
    integer n,i
    real randx(100000),theta,phi,omega
    real A1,A2,B1,B2,C1,C2,x,y,z,z1,z2,deg2rad

    deg2rad=atan2(1.0,1.0)/45.0
    theta=.0
    phi  =.0
    omega=.0
    A1=-sin(31.717*deg2rad)*cos(69.09*deg2rad);B1=-cos(31.717*deg2rad)*cos(69.09*deg2rad);C1=sin(31.717*deg2rad)*sin(69.09*deg2rad)
    A2=-sin(31.717*deg2rad)*cos(69.09*deg2rad);B2= cos(31.717*deg2rad)*cos(69.09*deg2rad);C2=sin(31.717*deg2rad)*sin(69.09*deg2rad)

    call randomreal(randx)
    n=0
    do i=first,last
100     n=n+1
        theta=69.0+randx(n)*21.0
        n=n+1
        phi=31.717*(randx(n)-0.5)*2.0
        x=sin(theta*deg2rad)*cos(phi*deg2rad)
        y=sin(theta*deg2rad)*sin(phi*deg2rad)
        z=cos(theta*deg2rad)
        z1=-(A1*x+B1*y)/C1
        z2=-(A2*x+B2*y)/C2
        if((z.ge..0).and.(z.le.z1).and.(z.le.z2)) then
            n=n+1
            omega=360.0*randx(n)
        else
            goto 100
        end if
        particle(i)%theta=theta
        particle(i)%phi  =phi
        particle(i)%omega=omega
    end do

end



