subroutine reconstruct_FSC
use common_icosproc
    implicit none
    integer hx,ky,lz,n1,n2,n3,i3,ir,tnfsize,mapsize,tmp
    real FSC1(0:maxR),FSC2(0:maxR),FSC(0:maxR),tmp2
    real,allocatable::a1(:),b1(:),a2(:),b2(:),wg1(:),wg2(:),fa(:),fb(:)

    open(11,file='oddtnf', form='unformatted',access='stream',status='old')
    open(12,file='evetnf',form='unformatted',access='stream',status='old')
    open(13,file='FSC.txt')
    print*,'calculating FSC...'
    tnfsize=(maxR+1)*(2*maxR+1)**2
    mapsize=FFTsize**3
    allocate(a1(tnfsize),b1(tnfsize),a2(tnfsize),b2(tnfsize),wg1(tnfsize),wg2(tnfsize))

    read(11) tmp,tmp,tmp2,tmp2,tmp
    read(11) a1
    read(11) b1
    read(11) wg1
    read(12) tmp,tmp,tmp2,tmp2,tmp
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
    a1=a1/wg1
    b1=b1/wg1
    deallocate(wg1)

    allocate(fa(mapsize),fb(mapsize))
    fa=0
    fb=0
    do lz=0,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                ir=int(sqrt(float(hx**2+ky**2+lz**2)))
                if(ir.le.maxR) then
                    i3=lz*(2*maxR+1)**2+(ky+maxR)*(2*maxR+1)+(hx+maxR)+1
                    n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
                    n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1  !n1-2*lz*FFTsize**2
                    if(mod(hx+ky+lz,2).eq.0) then
                        tmp=1
                    else
                        tmp=-1
                    end if
                    fa(n1)=a1(i3)*tmp
                    fb(n1)=b1(i3)*tmp
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                end if
            end do
        end do
    end do
    deallocate(a1,b1)

    print*,'perform 3D IFFT...'
    n1=FFTsize;n2=FFTsize;n3=FFTsize
    call FFT1d(fa,fb,n1*n2*n3,n1,n1,-1)
    print*,'IFFT 2'
    call FFT1d(fa,fb,n1*n2*n3,n2,n1*n2,-1)
    print*,'IFFT 3'
    call FFT1d(fa,fb,n1*n2*n3,n3,n1*n2*n3,-1)
    deallocate(fb)

    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        !tmp=(-1)**(hx+ky+lz)
        if(mod(hx+ky+lz,2).ne.0) fa(i3)=-fa(i3)
    end do
    call maphead
    open(21,file=trim(result3d),form='unformatted',access='stream')
    write(21) mrc

    if((proc3d%realspaceavg.eq.'Y').or.(proc3d%realspaceavg.eq.'y')) then
        call img3d_icosavg(fa,FFTsize,imgmask)
    end if
    print*,'writing 3D map...'
    do lz=1,FFTsize
        write(21) fa((lz-1)*FFTsize**2+1:lz*FFTsize**2)
    end do
    close(21)
    deallocate(fa)
    return
end



