
subroutine reconstruct_substep
use common_local_reconstruct
    implicit none
    integer i,i2,i3,j,hx,ky,lz,imgID0,no0,no1,imgparticle0,tmp1,tmp2,imgparticleID0
    integer usedparticle,signz,mapsize,tmaxR
    integer isym,lx(3),cubID(8)
    real theta,phi,omega,dr,rot(3,3),x0(3),x(3),x00(3)
    real ax(3),disx(8)
    real,allocatable::a(:),b(:),shift2d(:),img2d(:)
    real,allocatable::fa(:),fb(:),wg(:),img2dstck(:)

    open(21,file=trim(result3d), form='unformatted',access='stream',status='replace')
 
    tmaxR=2*maxR+1 
    mapsize=(tmaxR)**2*(maxR+1) 
    allocate(img2d(FFTsize0**2))
    allocate(a(FFTsize**2),b(FFTsize**2),shift2d(FFTsize**2))
    allocate(fa(mapsize),fb(mapsize),wg(mapsize))
    fa=.0
    fb=.0
    wg=.1
    imgID0=-100
    imgparticleID0=-100

    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        shift2d(i2)=(-1)**(hx+ky)
    end do

    no0=-1
    stckfile='abcdef'
    allocate(img2dstck(FFTsize0**2))
    if(proc2d%imgstck.eq.'y') then
        open(11,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(11) mrc
        do i=1,first-1
            call fseek(11,4*FFTsize0**2,1)
        end do
    end if

    do i=first,last
       ! print*,'i=',i
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        if(proc2d%imgstck.eq.'y') then
            read(11) img2d
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
            if(imgparticleID0.ne.particle(i)%imgparticleID.or.particle(i)%imgID.ne.imgID0) then
                img2d(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize0**2+1:particle(i)%imgparticleID*FFTsize0**2)
                call img2d_norm(img2d,FFTsize0)
                imgparticleID0=particle(i)%imgparticleID
            end if
          
        end if
        if(particle(i)%quality.eq.'g') then
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
               if((proc3d%Zhigh_df.eq.'y').or.(proc3d%Zhigh_df.eq.'Y'))then
                if(trim(sign).eq.'-')then
                 particle(i)%df1=particle(i)%df1-particle(i)%local_z
                 particle(i)%df2=particle(i)%df2-particle(i)%local_z
                 elseif(trim(sign).eq.'+')then
                 particle(i)%df1=particle(i)%df1+particle(i)%local_z
                 particle(i)%df2=particle(i)%df2+particle(i)%local_z
                endif
                    call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
                                apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                    else
                if(particle(i)%imgID.ne.imgID0) then
                    call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
                                apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                    imgID0=particle(i)%imgID
                end if
                end if
            end if

            call img2d_local_box(img2d,FFTsize0,particle(i)%local_hx0,particle(i)%local_ky0,a,FFTsize)
            call img2d_mask(a,FFTsize,imgmask,particle(i)%local_x,particle(i)%local_y)
            call img2d_norm(a,FFTsize)
            a=a*shift2d
            b=.0
            call img2d_FFT(a,b,FFTsize,FFTsize,1)
            call img2d_FT_centshift(a,b,FFTsize,particle(i)%local_x,particle(i)%local_y)
            theta=particle(i)%local_theta
            phi  =particle(i)%local_phi
            omega=particle(i)%local_omega
            call euler_mat(theta,phi,omega,rot)

            do ky=-maxR,maxR
                do hx=0,maxR
                    dr=sqrt(float(hx**2+ky**2))
                    if(dr.le.maxR-1) then
                        i2=(ky+FFTsize/2)*FFTsize+(hx+FFTsize/2)+1
                        x(1)=hx
                        x(2)=ky
                        x(3)=.0
                        x0=matmul(rot,x)
!************************************************************************************
                        x00=x0
                        do isym=1,Nsym
                            x=x00
                            x0=matmul(matsym(:,:,isym),x)
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
    write(21) FFTsize,maxR,apix,imgmask,tmp1,tmp2
    write(21) fa
    write(21) fb
    write(21) wg
    close(21)
    deallocate(a,b,fa,fb,wg,shift2d)
    return
end

