
subroutine reconstruct_C1
use common_icosproc
use FFTW3SPACE
    implicit none
    integer i,i2,j,hx,ky,lz,imgID0,no0,no1,n1,n2,imgparticle0
    integer usedparticle,signz,mapsize
    integer isym,lx(3),cubID(8)
    real theta,phi,omega,dr,rot(3,3),x0(3),x(3),x00(3),dcent
    real matsymc12(3,3,12),matsymc5a(3,3,5),matsymc6(3,3,6),matsymc3(3,3,3),matsymc15(3,3,15),matsymc4(3,3,4)
    real matsym(3,3,1),ax(3),disx(8),matsymc8(3,3,8),matsymc2(3,3,2),matsymc7(3,3,7)

    real,allocatable::a(:),b(:),shift2d(:)
    real,allocatable::fa(:),fb(:),wg(:),img2dstck(:)
    complex,allocatable::FFT3d(:)

    print*,'reconstruction'
    open(21,file=trim(result3d), form='unformatted',access='stream',status='replace')
    mapsize=FFTsize**3
    allocate(a(FFTsize**2),b(FFTsize**2),shift2d(FFTsize**2))
    allocate(fa(mapsize))
    allocate(fb(mapsize))
    allocate(wg(mapsize))
    !boundX=FFTsize/4.0
    fa=.0
    fb=.0
    wg=.01
    imgID0=-100
    imgparticle0=-100
    usedparticle=0
    matsym=.0
    matsym(1,1,1)=1;matsym(2,2,1)=1;matsym(3,3,1)=1
    matsymc12=.0
    matsymc12(1,1,1)=1;matsymc12(2,2,1)=1;matsymc12(3,3,1)=1
    matsymc15=.0
    matsymc15(1,1,1)=1;matsymc15(2,2,1)=1;matsymc15(3,3,1)=1
    matsymc6=.0
    matsymc6(1,1,1)=1;matsymc6(2,2,1)=1;matsymc6(3,3,1)=1
    matsymc3=.0
    matsymc3(1,1,1)=1;matsymc3(2,2,1)=1;matsymc3(3,3,1)=1
    matsymc4=.0
    matsymc4(1,1,1)=1;matsymc4(2,2,1)=1;matsymc4(3,3,1)=1
    matsymc5a=.0
    matsymc5a(1,1,1)=1;matsymc5a(2,2,1)=1;matsymc5a(3,3,1)=1
    matsymc2=.0
    matsymc2(1,1,1)=1;matsymc2(2,2,1)=1;matsymc2(3,3,1)=1
    matsymc8=.0
    matsymc8(1,1,1)=1;matsymc8(2,2,1)=1;matsymc8(3,3,1)=1
    matsymc7=.0
    matsymc7(1,1,1)=1;matsymc7(2,2,1)=1;matsymc7(3,3,1)=1
    if(sym.eq.'c12')  call icos2f_sym_matc12(matsymc12)
    if(sym.eq.'c15')   call icos2f_sym_matc15(matsymc15)
    if(sym.eq.'c5a')   call icos2f_sym_matc5a(matsymc5a)
    if(sym.eq.'c6')   call icos2f_sym_matc6(matsymc6)
    if(sym.eq.'c3') call icos2f_sym_matc3(matsymc3)
    if(sym.eq.'c4') call icos2f_sym_matc4(matsymc4)
    if(sym.eq.'c8') call icos2f_sym_matc8(matsymc8)
    if(sym.eq.'c2') call icos2f_sym_matc2(matsymc2)
    if(sym.eq.'c7') call icos2f_sym_matc7(matsymc7)
    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        shift2d(i2)=(-1)**(hx+ky)
    end do

    no0=-1
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))
    if(proc2d%imgstck.eq.'y') then
        open(11,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(11) mrc
        do i=1,first-1
            call fseek(11,4*FFTsize**2,1)
        end do
    end if

    do i=first,last
!        print*,'i=',i
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        if(proc2d%imgstck.eq.'y') then
            read(11) a
        else
            if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
                deallocate(img2dstck)
                stckfile=trim(particle(i)%stckfile)
                open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
                read(11) mrc
                allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
                read(11) img2dstck
                close(11)
            end if
            a(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        end if
        dcent=sqrt(particle(i)%x**2+particle(i)%y**2)
        if((particle(i)%PR.le.PR_threshold).and.(abs(particle(i)%x).le.boundX ).and.(abs(particle(i)%y).le.boundX)) then
            usedparticle=usedparticle+1
            call img2d_norm(a,FFTsize)
!        if(proc2d%imgstck.eq.'y') then
!        call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
!              apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
!        else
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
                if(particle(i)%imgID.ne.imgID0) then
                    call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
                                apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                    imgID0=particle(i)%imgID
                end if
            end if
 !       end if

            if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a,FFTsize,imgmask,particle(i)%x,particle(i)%y)
            a=a*shift2d
            b=.0
            call img2d_FFT(a,b,FFTsize,FFTsize,1)
            call img2d_FT_centshift(a,b,FFTsize,particle(i)%x,particle(i)%y)

            theta=particle(i)%theta
            phi  =particle(i)%phi
            omega=particle(i)%omega
            call euler_mat(theta,phi,omega,rot)
            do ky=-maxR,maxR
                do hx=0,maxR
                    dr=sqrt(float(hx**2+ky**2))
                    if(dr.le.maxR) then
                        i2=(ky+FFTsize/2)*FFTsize+(hx+FFTsize/2)+1
                        x(1)=hx
                        x(2)=ky
                        x(3)=.0
                        x0=matmul(rot,x)
!************************************************************************************
                        x00=x0
                       do isym=1,Nsym
                            x=x00
                            if(sym.eq.'c12')x0=matmul(matsymc12(:,:,isym),x)
                            if(sym.eq.'c15')x0=matmul(matsymc15(:,:,isym),x)
                            if(sym.eq.'c6')x0=matmul(matsymc6(:,:,isym),x)
                            if(sym.eq.'c3')x0=matmul(matsymc3(:,:,isym),x)
                            if(sym.eq.'c5a')x0=matmul(matsymc5a(:,:,isym),x)
                            if(sym.eq.'c4')x0=matmul(matsymc4(:,:,isym),x)
                            if(sym.eq.'c2')x0=matmul(matsymc2(:,:,isym),x)
                            if(sym.eq.'c8')x0=matmul(matsymc8(:,:,isym),x)
                            if(sym.eq.'c7')x0=matmul(matsymc7(:,:,isym),x)
                            if(x0(3).lt..0) then
                                x0=-x0
                                signz=-1
                            else
                                signz=1
                            end if
                            x0=x0+FFTsize/2
!************************************************************************************
                            lx=int(x0)
                            ax=x0-lx
                            disx(1)=(1.0-ax(1))*(1.0-ax(2))*(1.0-ax(3))
                            disx(2)=(    ax(1))*(1.0-ax(2))*(1.0-ax(3))
                            disx(3)=(1.0-ax(1))*(    ax(2))*(1.0-ax(3))
                            disx(4)=(    ax(1))*(    ax(2))*(1.0-ax(3))
                            disx(5)=(1.0-ax(1))*(1.0-ax(2))*(    ax(3))
                            disx(6)=(    ax(1))*(1.0-ax(2))*(    ax(3))
                            disx(7)=(1.0-ax(1))*(    ax(2))*(    ax(3))
                            disx(8)=(    ax(1))*(    ax(2))*(    ax(3))
                            cubID(1)=lx(3)*FFTsize**2+lx(2)*FFTsize+lx(1)+1
                            cubID(2)=cubID(1)+1
                            cubID(3:4)=cubID(1:2)+FFTsize
                            cubID(5:8)=cubID(1:4)+FFTsize**2
                            do j=1,8
                                fa(cubID(j))=fa(cubID(j))+a(i2)*disx(j)*ctf2d(i2)
                                fb(cubID(j))=fb(cubID(j))+b(i2)*disx(j)*ctf2d(i2)*signz
                                wg(cubID(j))=wg(cubID(j))+      disx(j)*ctf2d(i2)**2
                            end do
                        end do
                    end if
                end do
            end do
        end if
    end do

    lz=0
    do ky=1,maxR
        do hx=-maxR,maxR
            n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
            n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1
            dr=sqrt(float(hx**2+ky**2+lz**2))
            if(dr.le.maxR) then
                wg(n1)=wg(n1)+wg(n2)
                fa(n1)=(fa(n1)+fa(n2))/wg(n1)
                fb(n1)=(fb(n1)-fb(n2))/wg(n1)
                if(mod(hx+ky+lz,2).eq.0) then
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                else
                    fa(n1)=-fa(n1)
                    fb(n1)=-fb(n1)
                    fa(n2)= fa(n1)
                    fb(n2)=-fb(n1)
                end if
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
        n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
        n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1
        dr=sqrt(float(hx**2+ky**2+lz**2))
        if(dr.le.maxR) then
            wg(n1)=wg(n1)+wg(n2)
            fa(n1)=(fa(n1)+fa(n2))/wg(n1)
            fb(n1)=(fb(n1)-fb(n2))/wg(n1)
            if(mod(hx+ky+lz,2).eq.0) then
                fa(n2)= fa(n1)
                fb(n2)=-fb(n1)
            else
                fa(n1)=-fa(n1)
                fb(n1)=-fb(n1)
                fa(n2)= fa(n1)
                fb(n2)=-fb(n1)
            end if
        else
            fa(n1)=.0
            fb(n1)=.0
            fa(n2)=.0
            fb(n2)=.0
        end if
    end do

    hx=0
    n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
    fa(n1)=fa(n1)/wg(n1)
    fb(n1)=.0


    do lz=1,maxR
        do ky=-maxR,maxR
            do hx=-maxR,maxR
                n1=( lz+FFTsize/2)*FFTsize**2+( ky+FFTsize/2)*FFTsize+( hx+FFTsize/2)+1
                n2=(-lz+FFTsize/2)*FFTsize**2+(-ky+FFTsize/2)*FFTsize+(-hx+FFTsize/2)+1
                dr=sqrt(float(hx**2+ky**2+lz**2))
                if(dr.le.maxR) then
                    fa(n1)=fa(n1)/wg(n1)
                    fb(n1)=fb(n1)/wg(n1)
                    if(mod(hx+ky+lz,2).eq.0) then
                        fa(n2)= fa(n1)
                        fb(n2)=-fb(n1)
                    else
                        fa(n1)=-fa(n1)
                        fb(n1)=-fb(n1)
                        fa(n2)= fa(n1)
                        fb(n2)=-fb(n1)
                    end if
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

    allocate(FFT3d(mapsize))
    FFT3d=cmplx(fa,-fb)
    deallocate(fb)
    call FFTNEW(4)
    print *,'Doing invert 3D FFT!'
    call FFTcmplx3d_B(FFT3d,FFT3d,FFTsize,FFTsize,FFTsize)
    call shift3dComplex(FFT3d,FFTsize,FFTsize,FFTsize)
    fa=real(FFT3d)
    deallocate(FFT3d)

    call maphead
    write(21) mrc

    if((proc3d%realspaceavg.eq.'Y').or.(proc3d%realspaceavg.eq.'y')) then
        call img3d_icosavg(fa,FFTsize,imgmask)
    end if
    write(21) fa
    print*,'usedparticle:',usedparticle

!    do lz=1,FFTsize
!        write(21) fa((lz-1)*FFTsize**2+1:lz*FFTsize**2)
!    end do
    if(proc2d%imgstck.eq.'y')  close(11)
    close(21)
    deallocate(fa,a,b,shift2d)
    return
end

