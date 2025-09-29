
subroutine guass_regression(z,dx,dy,cc)
    implicit none
    integer i,j
    real dx,dy,c10,c01,c11,c20,c02,c00,cc,z(-1:1,-1:1)

    c10=.0;c01=.0;c11=.0;c20=.0;c02=.0;c00=.0

    do j=-1,1
        do i=-1,1
            c10=c10+float(i)*log(z(i,j))
        end do
    end do
    c10=c10/6.0

    do j=-1,1
        do i=-1,1
            c01=c01+float(j)*log(z(i,j))
        end do
    end do
    c01=c01/6.0

    do j=-1,1
        do i=-1,1
            c11=c11+float(i*j)*log(z(i,j))
        end do
    end do
    c11=c11/4.0

    do j=-1,1
        do i=-1,1
            c20=c20+float(3*i**2-2)*log(z(i,j))
        end do
    end do
    c20=c20/6.0

    do j=-1,1
        do i=-1,1
            c02=c02+float(3*j**2-2)*log(z(i,j))
        end do
    end do
    c02=c02/6.0

    do j=-1,1
        do i=-1,1
            c00=c00+float(5-3*i**2-3*j**2)*log(z(i,j))
        end do
    end do
    c00=c00/9.0

    dx=(c11*c01-2.0*c10*c02)/(4.0*c20*c02-c11**2)
    dy=(c11*c10-2.0*c01*c20)/(4.0*c20*c02-c11**2)
    cc=exp(c00-c20*dx**2-c11*dx*dy-c02*dy**2)

    return
end
