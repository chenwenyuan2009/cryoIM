subroutine reconstruct_local_substep
use common_icosproc
use FFTW3SPACE
    implicit none
    integer i,i2,j,hx,ky,lz,imgID0,no0,no1,n1,n2,imgparticle0
    integer usedparticle,signz,mapsize,ix,iy,ip1,ip2,tmaxR
    integer isym,lx(3),cubID(8),x10,y10,tmp1,tmp2
    real theta,phi,omega,dr,rot(3,3),x0(3),x(3),x00(3),dcent
    real matsymc12(3,3,12),matsymc5a(3,3,5),matsymc6(3,3,6),matsymc3(3,3,3),matsymc15(3,3,15),matsymc4(3,3,4)
    real matsym(3,3,60),ax(3),disx(8),matsymc8(3,3,8),matsymc2(3,3,2)

    real,allocatable::a(:),b(:),shift2d(:)
    real,allocatable::fa(:),fb(:),wg(:),img2dstck(:),img2do(:)
    complex,allocatable::FFT3d(:)

    print*,'reconstruction'
    open(21,file=trim(result3d), form='unformatted',access='stream',status='replace')
    tmaxR=2*maxR+1
    mapsize=(tmaxR)**2*(maxR+1)
    allocate(a(loc_box**2),b(loc_box**2),shift2d(loc_box**2))
    allocate(fa(mapsize),img2do(FFTsize**2))
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
    if(sym.eq.'c12')  call icos2f_sym_matc12(matsymc12)
    if(sym.eq.'c15')   call icos2f_sym_matc15(matsymc15)
    if(sym.eq.'c5a')   call icos2f_sym_matc5a(matsymc5a)
    if(sym.eq.'c6')   call icos2f_sym_matc6(matsymc6)
    if(sym.eq.'c3') call icos2f_sym_matc3(matsymc3)
    if(sym.eq.'c4') call icos2f_sym_matc4(matsymc4)
    if(sym.eq.'c8') call icos2f_sym_matc8(matsymc8)
    if(sym.eq.'c2') call icos2f_sym_matc2(matsymc2)
    do i2=1,loc_box**2
        ky=(i2-1)/loc_box
        hx=mod(i2-1,loc_box)
        shift2d(i2)=(-1)**(hx+ky)
    end do

    no0=-1
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))
    do i=first,last
       if((proc2d%check.eq.'y').or.(proc2d%check.eq.'Y')) print*,i
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
            if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
                deallocate(img2dstck)
                stckfile=trim(particle(i)%stckfile)
                open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
                read(11) mrc
                allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
                read(11) img2dstck
                close(11)
            end if
            img2do(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
            x10=int(particle(i)%x0)+FFTsize/2
            y10=int(particle(i)%y0)+FFTsize/2            
            do iy=y10-loc_box/2,y10+loc_box/2-1
                do ix=x10-loc_box/2,x10+loc_box/2-1
                    ip1=iy*FFTsize+ix+1
                    ip2=(iy-y10+loc_box/2)*loc_box+ix-x10+loc_box/2+1
                    a(ip2) =img2do(ip1)
                end do
            end do

        !dcent=sqrt(particle(i)%x**2+particle(i)%y**2)
        if((particle(i)%PR.le.PR_threshold).and.(abs(particle(i)%x).le.boundX ).and.(abs(particle(i)%y).le.boundX)) then
            usedparticle=usedparticle+1
            call img2d_norm(a,loc_box)
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
                if(particle(i)%imgID.ne.imgID0) then
                    call getCTF2d(loc_box,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
                                apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                    imgID0=particle(i)%imgID
                end if
            end if

            if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a,loc_box,imgmask,particle(i)%x,particle(i)%y)
            a=a*shift2d
            b=.0
            call img2d_FFT(a,b,loc_box,loc_box,1)
            call img2d_FT_centshift(a,b,loc_box,particle(i)%x,particle(i)%y)

            theta=particle(i)%theta
            phi  =particle(i)%phi
            omega=particle(i)%omega
            call euler_mat(theta,phi,omega,rot)

            do ky=-maxR,maxR
                do hx=0,maxR
                    dr=sqrt(float(hx**2+ky**2))
                    if(dr.le.maxR-1) then
                        i2=(ky+loc_box/2)*loc_box+(hx+loc_box/2)+1
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
                            if(sym.eq.'c1')x0=matmul(matsym(:,:,isym),x)
                            if(x0(3).lt..0) then
                                x0=-x0
                                signz=-1
                            else
                                signz=1
                            end if
                            x0(1:2)=x0(1:2)+maxR
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
                            cubID(1)=lx(3)*tmaxR**2+lx(2)*tmaxR+lx(1)+1
                            cubID(2)=cubID(1)+1
                            cubID(3:4)=cubID(1:2)+tmaxR
                            cubID(5:8)=cubID(1:4)+tmaxR**2
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
    tmp1=0
    tmp2=0
    write(21) loc_box,maxR,apix,imgmask,tmp1,tmp2
    write(21) fa
    write(21) fb
    write(21) wg
    close(21)
    deallocate(a,b,shift2d,fa,img2do,fb,wg)

    return
end

     