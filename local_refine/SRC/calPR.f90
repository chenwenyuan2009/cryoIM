subroutine calPR(a1,a2,nx,minR,maxR,PR)
    implicit none
    integer nx,hx,ky,i2
    real minR,maxR,PR,PR1,PR2,dr,amp1,amp2,phase1,phase2,diff,pi,Tpi
    real a1(nx**2),a2(nx**2)
    real,allocatable::b1(:),b2(:)
    allocate(b1(nx**2),b2(nx**2))

    pi=4.0*atan2(1.0,1.0)
    Tpi=2.0*pi
    b1=.0;b2=.0
    call img2d_FFT(a1,b1,nx,nx,1)
    call img2d_FFT(a2,b2,nx,nx,1)
    PR1=.0;PR2=.0
    do i2=1,nx**2
        ky=(i2-1)/nx-nx/2
        hx=mod(i2-1,nx)-nx/2
        dr=sqrt(float(hx**2+ky**2))
        if((dr.ge.minR).and.(dr.le.maxR)) then
            amp1=sqrt(a1(i2)**2+b1(i2)**2)
            amp2=sqrt(a2(i2)**2+b2(i2)**2)
            phase1=atan2(b1(i2),a1(i2))
            phase2=atan2(b2(i2),a2(i2))
            diff=abs(phase2-phase1)
            diff=amod(diff,Tpi)
            if(diff.gt.pi) diff=Tpi-diff
            amp1=0.5*(amp1+amp2)
            PR1=diff*amp1+PR1
            PR2=amp1+PR2
        end if
    end do

    if(PR2.gt.1.0e-10) then
        PR=PR1/PR2*180.0/pi
    else
        PR=1000.
    end if

    return
    end

        