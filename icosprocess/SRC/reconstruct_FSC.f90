subroutine reconstruct_FSC
use FFTW3SPACE
use common_icosproc
    implicit none
    integer hx,ky,lz,n1,n2,i3,i3a,i3b,ir,tnfsize,mapsize,tmp
    real tmp2
    real,allocatable::FSC1(:),FSC2(:),FSC(:),a1(:),b1(:),a2(:),b2(:),wg1(:),wg2(:),fa(:),fb(:)
    complex,allocatable::FFT3d(:)

    allocate(FSC1(0:maxR),FSC2(0:maxR),FSC(0:maxR))

    open(11,file='oddtnf', form='unformatted',access='stream',status='old')
    open(12,file='evetnf',form='unformatted',access='stream',status='old')
    open(13,file='FSC.txt')
    print*,'calculating FSC...'
    tnfsize=(maxR+1)*(2*maxR+1)**2
    mapsize=FFTsize**3
    allocate(a1(tnfsize),b1(tnfsize),a2(tnfsize),b2(tnfsize),wg1(tnfsize),wg2(tnfsize))

    read(11) tmp,tmp,tmp2,tmp2,tmp,tmp
    read(11) a1
    read(11) b1
    read(11) wg1
    read(12) tmp,tmp,tmp2,tmp2,tmp,tmp
    read(12) a2
    read(12) b2
    read(12) wg2
    close(11)
    close(12)
    wg1=wg1+1.0
    wg2=wg2+1.0
    a1=a1/wg1
    b1=b1/wg1
    a2=a2/wg2
    b2=b2/wg2

    FSC=.0
    FSC1=1.0E-10
    FSC2=1.0E-10
    do lz=0,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                ir=int(sqrt(float(hx**2+ky**2+lz**2)))
                if(ir.le.maxR) then
                    i3=lz*(2*maxR+1)**2+(ky+maxR)*(2*maxR+1)+(hx+maxR)+1
                    FSC(ir) =FSC(ir)+(a1(i3)*a2(i3)+b1(i3)*b2(i3))
                    FSC1(ir)=FSC1(ir)+(a1(i3)**2+b1(i3)**2)
                    FSC2(ir)=FSC2(ir)+(a2(i3)**2+b2(i3)**2)
                end if
            end do
        end do
    end do

    do ir=0,maxR
        FSC(ir)=FSC(ir)/sqrt(FSC1(ir)*FSC2(ir))
        write(13,*) ir,FSC(ir),sqrt(2.0*FSC(ir)/(1.0+FSC(ir)))
        if(FSC(ir).le..0) then
            FSC(ir)=.0
        else
            FSC(ir)=sqrt(2.0*FSC(ir)/(1.0+FSC(ir)))
        end if
    end do
    close(13)

    a1=a1*wg1
    b1=b1*wg1
    a2=a2*wg2
    b2=b2*wg2
    a1=a1+a2
    b1=b1+b2
    wg1=wg1+wg2-2.0
    deallocate(a2,b2,wg2)
    do lz=0,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                ir=int(sqrt(float(hx**2+ky**2+lz**2)))
                if(ir.le.maxR) then
                    i3=lz*(2*maxR+1)**2+(ky+maxR)*(2*maxR+1)+(hx+maxR)+1
                    a1(i3) =a1(i3) *FSC(ir)
                    b1(i3) =b1(i3) *FSC(ir)
                    wg1(i3)=wg1(i3)*FSC(ir)
                else
                    a1(i3)=.0
                    b1(i3)=.0
                end if
            end do
        end do
    end do
    wg1=wg1+1.0
!    a1=a1/wg1
!    b1=b1/wg1
!    deallocate(wg1)


    allocate(fa(mapsize),fb(mapsize))
    fa=0
    fb=0
    lz=0
    do ky=1,maxR
        do hx=-maxR,maxR
            ir=int(sqrt(float(hx**2+ky**2+lz**2)))
            if(ir.le.maxR) then
                i3a= lz*(2*maxR+1)**2+( ky+maxR)*(2*maxR+1)+( hx+maxR)+1
                i3b=-lz*(2*maxR+1)**2+(-ky+maxR)*(2*maxR+1)+(-hx+maxR)+1
                n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
                n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1
                wg1(i3a)=wg1(i3a)+wg1(i3b)
                fa(n1)=(a1(i3a)+a1(i3b))/wg1(i3a)
                fb(n1)=(b1(i3a)-b1(i3b))/wg1(i3a)
                if(mod(hx+ky+lz,2).eq.0) then
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                else
                    fa(n1)=-fa(n1)
                    fb(n1)=-fb(n1)
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                end if
            end if
        end do
    end do

    ky=0
    do hx=1,maxR
        ir=int(sqrt(float(hx**2+ky**2+lz**2)))
        if(ir.le.maxR) then
            i3a= lz*(2*maxR+1)**2+( ky+maxR)*(2*maxR+1)+( hx+maxR)+1
            i3b=-lz*(2*maxR+1)**2+(-ky+maxR)*(2*maxR+1)+(-hx+maxR)+1
            n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
            n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1
            wg1(i3a)=wg1(i3a)+wg1(i3b)
            fa(n1)=(a1(i3a)+a1(i3b))/wg1(i3a)
            fb(n1)=(b1(i3a)-b1(i3b))/wg1(i3a)
            if(mod(hx+ky+lz,2).eq.0) then
                fa(n2)= fa(n1)
                fb(n2)=-fb(n1)
            else
                fa(n1)=-fa(n1)
                fb(n1)=-fb(n1)
                fa(n2)= fa(n1)
                fb(n2)=-fb(n1)
            end if
        end if
    end do

    hx=0
    i3a= lz*(2*maxR+1)**2+( ky+maxR)*(2*maxR+1)+( hx+maxR)+1
    n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
    fa(n1)=a1(i3a)/wg1(i3a)
    fb(n1)=.0

    do lz=1,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                ir=int(sqrt(float(hx**2+ky**2+lz**2)))
                if(ir.le.maxR) then
                    i3=lz*(2*maxR+1)**2+(ky+maxR)*(2*maxR+1)+(hx+maxR)+1
                    n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
                    n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1  !n1-2*lz*FFTsize**2
                    tmp=-1
                    if(mod(hx+ky+lz,2).eq.0) tmp=1
                    fa(n1)=a1(i3)/wg1(i3)*tmp
                    fb(n1)=b1(i3)/wg1(i3)*tmp
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                end if
            end do
        end do
    end do
    deallocate(a1,b1)

    allocate(FFT3d(mapsize))
    FFT3d=cmplx(fa,fb)
    deallocate(fb)
    call FFTNEW(4)
    print *,'Doing invert 3D FFT!'
    call FFTcmplx3d_B(FFT3d,FFT3d,FFTsize,FFTsize,FFTsize)
    call shift3dComplex(FFT3d,FFTsize,FFTsize,FFTsize)
    fa=real(FFT3d)
    deallocate(FFT3d)

    open(21,file=trim(result3d),form='unformatted',access='stream')
    write(21) mrc

    if((proc3d%realspaceavg.eq.'Y').or.(proc3d%realspaceavg.eq.'y')) then
        call img3d_icosavg(fa,FFTsize,imgmask)
    end if
    write(21) fa

!    print*,'writing 3D map...'
!    do lz=1,FFTsize
!        write(21) fa((lz-1)*FFTsize**2+1:lz*FFTsize**2)
!    end do
    close(21)
    deallocate(fa)
    return
end



