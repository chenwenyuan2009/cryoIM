subroutine algin_findvertex
use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1
    integer ix0,iy0,ix,iy,n1,ip1,ip2,x0,y0
    real dx,dy,cc,dn,dm,RotMat(3,3)
    real ort60(3,60),cent60(3,60),Cmax,z(-1:1,-1:1),initiald(3,1),newd(3,1)
    real,allocatable::a0(:),a1(:),a2(:),a(:),a3(:)
    real,allocatable::core_FTr(:),core_FTi(:)
    real,allocatable::loc_img2d(:,:),img2dstck(:),loc_proj2d(:,:)
    real(kind=8) Pr,phaseRedisual
    allocate(a0(FFTsize**2),a1(FFTsize**2),a2(loc_box**2),a3(loc_box**2))
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),a(FFTsize**2),img2dstck(FFTsize**2))
    allocate(loc_img2d(loc_box**2,60),loc_proj2d(loc_box**2,60))

    print*,'reading 3D map...'
    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    if(proc2d%imgstck.eq.'y')then
            open(12,file=trim(imgstck),form='unformatted',access='stream',status='old')
            read(12) mrc
    endif
    open(21,file=trim(newortfile))
    read(11) mrc
    do i=1,FFTsize
      read(11) a
      do j=1,FFTsize**2
        core_FTr((i-1)*FFTsize**2+j)=a(j)
      end do
    end do

    close(11)

    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            core_FTr(i3)  =-core_FTr(i3)
        end if
    end do
    core_FTi  =.0
    stckfile='abcd'
    dn=0.52573*particle_radius/apix
    dm=0.85065*particle_radius/apix
    initiald(1,1)=0.0;initiald(2,1)=dm;initiald(3,1)=dn

    print*,'preform 3D FFT...'
    call img3d_FFT(core_FTr,core_FTi,FFTsize,FFTsize,FFTsize,1)
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            core_FTr(i3)=-core_FTr(i3)
            core_FTi(i3)=-core_FTi(i3)
        end if
    end do

    no0=-1
    imgID0=-100

    do i=1,first-1
        call fseek(12,4*FFTsize**2,1)
    end do
    do i=first,last
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
       if(proc2d%imgstck.eq.'y')then
          read(12)a0
       else
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(12,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(12) mrc
            allocate(img2dstck(FFTsize**2*mrc%nz))
            read(12)img2dstck
            close(12)
        end if
        a0(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        endif
        call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)
        call img2d_lowpass(a0,FFTsize,maxR)
        !call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        call img2d_norm(a0,FFTsize)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if
        call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,particle(i)%theta,particle(i)%phi,particle(i)%omega,a1)
        call img2d_applyCTF(a1,FFTsize,ctf2d)
        call img2d_lowpass(a1,FFTsize,maxR)
        !call img2d_mask(a1,FFTsize,imgmask,.0,.0)
        call icos_equalort(particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60)
        !print*,particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60(:,1)
        do j=1,60

            ort60(2,j)=ort60(2,j)
            ort60(3,j)=ort60(3,j)
            call euler_mat_vertex(ort60(1,j),ort60(2,j),ort60(3,j),RotMat)
             newd=matmul(RotMat,initiald)

            x0=int(newd(1,1))+FFTsize/2
            y0=int(newd(2,1))+FFTsize/2

            do iy=y0-loc_box/2,y0+loc_box/2-1
                do ix=x0-loc_box/2,x0+loc_box/2-1
                    ip1=iy*FFTsize+ix+1
                    ip2=(iy-y0+loc_box/2)*loc_box+ix-x0+loc_box/2+1
                    loc_img2d(ip2,j) =a0(ip1)
                    loc_proj2d(ip2,j)=a1(ip1)
                end do
            end do

            a2=loc_img2d(:,j)
            a3=loc_proj2d(:,j)
            call img2d_mask(a2,loc_box,imgmask,.0,.0)
            call img2d_mask(a3,loc_box,imgmask,.0,.0)
            call img2d_norm(a2,loc_box)
            call img2d_norm(a3,loc_box)
!            call cross_correlation(a2,a3,loc_box)
            Pr=phaseRedisual(a2,a3,loc_box)
            if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) then
                Nmax=maxloc(a2)-1
                iy0=Nmax(1)/FFTsize
                ix0=mod(Nmax(1),FFTsize)
                do ky=-1,1
                    do hx=-1,1
                        ix=ix0+hx
                        iy=iy0+ky
                        n1=iy*FFTsize+ix+1
                        z(hx,ky)=a2(n1)
                    end do
                end do
                call guass_regression(z,dx,dy,cc)
                cent60(1,j)=float(ix0-FFTsize/2)+dx
                cent60(2,j)=float(iy0-FFTsize/2)+dy
                cent60(3,j)=cc
            else
                cent60(1,j)=particle(i)%x
                cent60(2,j)=particle(i)%y
                !cent60(3,j)=a2(loc_box**2/2+loc_box/2+1)
                cent60(3,j)=Pr
            end if

        end do
        Nmax=minloc(cent60(3,:))
        if(cent60(3,Nmax(1)).ge.1) cent60(3,Nmax(1))=1
        particle(i)%theta=ort60( 1,Nmax(1))
        particle(i)%phi  =ort60( 2,Nmax(1))
        particle(i)%omega=ort60( 3,Nmax(1))
        particle(i)%x    =cent60(1,Nmax(1))
        particle(i)%y    =cent60(2,Nmax(1))
        particle(i)%PR   =cent60(3,Nmax(1))

        print*,i,Nmax,cent60(:,Nmax(1))

        Cmax=acos(Cmax)*45.0/atan2(1.0,1.0)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,   particle(i)%astigang,Nmax

    end do
    deallocate(a,a0,a1,a2,a3,loc_img2d,loc_proj2d)
    deallocate(core_FTr,core_FTi)
    close(21)
    if(proc2d%imgstck.eq.'y') close(12)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5)
    return
end

