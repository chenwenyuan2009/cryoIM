
!相对于reconstruct2减小傅立叶空间a,b数组的大小


subroutine merge_reconstruct_test
use common_merge_icosprocess_subrec
use FFTW3SPACE
    implicit none
    integer i,hx,ky,lz,n1,n2,n3,n,tmaxR,no0,no1
    integer tmp,dr
    real,allocatable::fa(:),fb(:),wg(:),fa1(:),fb1(:),wg1(:),fa0(:),fb0(:),wg0(:)
    complex,allocatable::FFT3d(:)

    open(21,file=trim(par(nfile)), form='unformatted',access='stream',status='replace')

    tmaxR=2*maxR+1
    allocate(fa0(tmaxR**2*(maxR+1)),fb0(tmaxR**2*(maxR+1)),wg0(tmaxR**2*(maxR+1)))
    allocate(fa(tmaxR**3),fb(tmaxR**3),wg(tmaxR**3),fa1(tmaxR**3),fb1(tmaxR**3),wg1(tmaxR**3))
    
    fa=.0
    fb=.0
    wg=0.01
    do i=1,nfile-1
        open(11,file=trim(par(i)),form='unformatted',access='stream',status='old')
        print*,'reading data form ',trim(par(i))
        read(11) FFTsize,maxR,apix,imgmask,nFSC,nicosavg
        read(11) fa0
        read(11) fb0
        read(11) wg0
        fa(tmaxR**2*maxR+1:tmaxR**3)=fa(tmaxR**2*maxR+1:tmaxR**3)+fa0(:)
        fb(tmaxR**2*maxR+1:tmaxR**3)=fb(tmaxR**2*maxR+1:tmaxR**3)+fb0(:)
        wg(tmaxR**2*maxR+1:tmaxR**3)=wg(tmaxR**2*maxR+1:tmaxR**3)+wg0(:)
        close(11)
    end do
    deallocate(wg0,fa0,fb0,wg1,fa1,fb1)

    lz=0
    do ky=1,maxR
        do hx=-maxR,maxR
            n1=( lz+maxR)*tmaxR**2+( ky+maxR)*tmaxR+( hx+maxR)+1
            n2=(-lz+maxR)*tmaxR**2+(-ky+maxR)*tmaxR+(-hx+maxR)+1
            dr=sqrt(float(hx**2+ky**2+lz**2))
            if(dr.le.maxR) then
                wg(n1)=wg(n1)+wg(n2)
                tmp=-1
                if(mod(hx+ky+lz,2).eq.0) tmp=1
                fa(n1)=(fa(n1)+fa(n2))/wg(n1)*tmp
                fb(n1)=(fb(n1)-fb(n2))/wg(n1)*tmp
                fa(n2)= fa(n1)
                fb(n2)=-fb(n1)
            else
                fa(n1)=.0
                fb(n1)=.0
                fa(n2)=.0
                fb(n2)=.0
            end if
        end do
    end do

    ky=0
    do hx=1,maxR
        n1=( lz+maxR)*tmaxR**2+( ky+maxR)*tmaxR+( hx+maxR)+1
        n2=(-lz+maxR)*tmaxR**2+(-ky+maxR)*tmaxR+(-hx+maxR)+1
        dr=sqrt(float(hx**2+ky**2+lz**2))
        if(dr.le.maxR) then
            wg(n1)=wg(n1)+wg(n2)
            tmp=-1
            if(mod(hx+ky+lz,2).eq.0) tmp=1
            fa(n1)=(fa(n1)+fa(n2))/wg(n1)*tmp
            fb(n1)=(fb(n1)-fb(n2))/wg(n1)*tmp
            fa(n2)= fa(n1)
            fb(n2)=-fb(n1)
        else
            fa(n1)=.0
            fb(n1)=.0
            fa(n2)=.0
            fb(n2)=.0
        end if
    end do

    hx=0
    n1=( lz+maxR)*tmaxR**2+( ky+maxR)*tmaxR+( hx+maxR)+1
    fa(n1)=fa(n1)/wg(n1)
    fb(n1)=.0

    do lz=1,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                n1=( lz+maxR)*tmaxR**2+( ky+maxR)*tmaxR+( hx+maxR)+1
                n2=(-lz+maxR)*tmaxR**2+(-ky+maxR)*tmaxR+(-hx+maxR)+1
                dr=sqrt(float(hx**2+ky**2+lz**2))
                if(dr.le.maxR) then
                    tmp=-1
                    if(mod(hx+ky+lz,2).eq.0) tmp=1
                    fa(n1)=fa(n1)*tmp/wg(n1)
                    fb(n1)=fb(n1)*tmp/wg(n1)
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                else
                    fa(n1)=.0
                    fb(n1)=.0
                    fa(n2)=.0
                    fb(n2)=.0
                end if
            end do
        end do
    end do
    deallocate(wg)

    allocate(fa1(FFTsize**3))
    fa1=.0
    do i=1,tmaxR**3
        lz=(i-1)/tmaxR**2-maxR+FFTsize/2
        ky=mod(i-1,tmaxR**2)/tmaxR-maxR+FFTsize/2
        hx=mod(mod(i-1,tmaxR**2),tmaxR)-maxR+FFTsize/2
        n=lz*FFTsize**2+ky*FFTsize+hx+1
        fa1(n)=fa(i)
    end do
    deallocate(fa)

    allocate(fb1(FFTsize**3))
    fb1=.0
    do i=1,tmaxR**3
        lz=(i-1)/tmaxR**2-maxR+FFTsize/2
        ky=mod(i-1,tmaxR**2)/tmaxR-maxR+FFTsize/2
        hx=mod(mod(i-1,tmaxR**2),tmaxR)-maxR+FFTsize/2
        n=lz*FFTsize**2+ky*FFTsize+hx+1
        fb1(n)=fb(i)
    end do
    deallocate(fb)

    allocate(FFT3d(FFTsize**3))
    FFT3d=cmplx(fa1,-fb1)
    deallocate(fb1)
    call FFTNEW(4)
    print *,'Doing invert 3D FFT!'
    call FFTcmplx3d_B(FFT3d,FFT3d,FFTsize,FFTsize,FFTsize)
    call shift3dComplex(FFT3d,FFTsize,FFTsize,FFTsize)
    fa1=real(FFT3d)
    deallocate(FFT3d) 

    if(nicosavg.eq.1) call img3d_icosavg(fa1,FFTsize,imgmask)
    call maphead
    write(21) mrc
    write(21) fa1
    close(21)

    deallocate(fa1)
    return
end


