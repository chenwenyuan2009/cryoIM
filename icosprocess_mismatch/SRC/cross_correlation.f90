subroutine cross_correlation(a1,a2,nx)
    implicit none
    integer nx,i,hx,ky
    real x,y
    real a1(nx**2),a2(nx**2)
    real b1(nx**2),b2(nx**2)

    b1=.0
    b2=.0
    call img2d_FFT(a1,b1,nx,nx,1)
    call img2d_FFT(a2,b2,nx,nx,1)

    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        x=a1(i)*a2(i)+b1(i)*b2(i)
        y=b1(i)*a2(i)-a1(i)*b2(i)
        a1(i)=x*(-1)**(hx+ky)
        b1(i)=y*(-1)**(hx+ky)
    end do

    call img2d_FFT(a1,b1,nx,nx,-1)
    return
end

function phaseRedisual(a1,a2,locaSize) result(pr)
    use common_icosproc
    implicit none
    integer i,n1,n2,h,k,locaSize,maxR1
    complex f1,f2
    real x,y,Amp1,amp2,phase1,phase2,diff
    real R,pi
    real(kind=8) Pr,ampweight,mean1,mean2
    real a1(locaSize*locaSize),a2(locaSize*locaSize)
    real a11(Padsize*Padsize)
    real a12(Padsize*Padsize)
    real b1(Padsize**2),b2(Padsize**2)
    pi=3.1415926
    mean1=0.0
    mean2=0.0

    do i=1,locaSize**2
        mean1=mean1+a1(i)
        mean2=mean2+a2(i)
    end do

    mean1=mean1/float(locaSize**2)
    mean2=mean2/float(locaSize**2)

    a11=mean1
    a12=mean2

    IF(Padsize.NE.locaSize) THEN
     call PADFFT(a1,a11,locaSize,Padsize)
     call PADFFT(a2,a12,locaSize,Padsize)
    else
     a11=a1
     a12=a2
    END IF

!    open(211,file='padtest.mrc',form='unformatted',access='stream')
!    mrc%nz=2
!    mrc%nx=Padsize
!    mrc%ny=Padsize
!    write(211) mrc
!    write(211) a11
!    write(211) a12
!    close(211)

     do i=1,(Padsize)**2
        k=(i-1)/(Padsize)
        h=mod(i-1,Padsize)
        R=sqrt(k**2.0+h**2.0)
        a12(i)=a12(i)*(-1)**(h+k)
        a11(i)=a11(i)*(-1)**(h+k)
    end do

    n1=Padsize
    n2=Padsize
    b1=.0
    b2=.0

    call FFT1d(a11,b1,n1*n2,n1,n1,1)
    call FFT1d(a11,b1,n1*n2,n2,n1*n2,1)
    call FFT1d(a12,b2,n1*n2,n1,n1,1)
    call FFT1d(a12,b2,n1*n2,n2,n1*n2,1)

    Pr=0.0
    ampweight=0.0
    maxR1=maxR*Padsize/locaSize
    do i=1,Padsize**2
        k=(i-1)/Padsize-Padsize/2
        h=mod(i-1,Padsize)-Padsize/2
        R=sqrt(k**2.0+h**2.0)
        if(R.gt. minR .and. R .lt.maxR1) then
          Amp1=sqrt(a11(i)**2+b1(i)**2)
          phase1=atan2(b1(i),a11(i))
          Amp2=sqrt(a12(i)**2+b2(i)**2)
          phase2=atan2(b2(i),a12(i))
          diff=abs(phase1-phase2)
          if(diff.gt.Pi) diff=2*PI-diff
          pr=pr+diff*(amp1+amp2)
          ampweight=ampweight+(amp1+amp2)
        end if
    end do

     pr=1.0-pr/ampweight/pi

end function


subroutine PADFFT(a1,a2,localsize,padsize1)
   use common_icosproc
   implicit none
   integer localsize,padsize1,i,j,ky,hx
   real a1(localsize**2)
   real a2(padsize1**2)
   real R

   do i=1,padsize1**2
      ky=(i-1)/Padsize1-Padsize1/2
      hx=mod(i-1,Padsize1)-Padsize1/2
      R=sqrt(ky**2.0+hx**2.0)
      if(R.lt.localsize/2) then
        ky=ky+localsize/2
        hx=hx+localsize/2
        j=ky*localsize+hx
        a2(i)=a1(j)
      end if
   end do

end subroutine


