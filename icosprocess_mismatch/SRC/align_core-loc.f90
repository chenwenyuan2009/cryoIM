subroutine algin_core_loc
use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1,NNmax(1),nproc
    integer ix0,iy0,ix,iy,n1,ip1,ip2,x0,y0,ix00,iy00,jj,ii
    real dx,dy,cc,dn,dm,RotMat(3,3),delts,dang,dtran,ortcent(5,11),dfcent(3,7),dfpr(7)
    real ort60(3),cent60(3),Cmax,z(-1:1,-1:1),initiald(3,1),newd(3,1),cent11(3,11)
    real dforigin(3)
    real,allocatable::a0(:),a1(:),a2(:),a(:),a3(:),a4(:)
    real,allocatable::core_FTr(:),core_FTi(:),img2dstck(:),loc_img2d2(:),loc_proj2d2(:)
    real,allocatable::loc_img2d(:,:),loc_proj2d(:,:),loc_img2d1(:,:),loc_proj2d1(:,:)
    real(kind=8) Pr,phaseRedisual

    allocate(a0(FFTsize**2),a1(FFTsize**2),a2(loc_box**2),a3(loc_box**2),loc_img2d2(loc_box**2),loc_proj2d2(loc_box**2))
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),a(FFTsize**2),img2dstck(FFTsize**2),a4(FFTsize**2))
    allocate(loc_img2d(loc_box**2,5),loc_proj2d(loc_box**2,5),loc_img2d1(loc_box**2,11),loc_proj2d1(loc_box**2,11))

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
    cent60=.0
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            core_FTr(i3)  =-core_FTr(i3)
        end if
    end do
    core_FTi  =.0
    delts=.0

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
    stckfile='abcd'
    if(proc2d%imgstck.eq.'y')then
    do i=1,first-1
        call fseek(12,4*FFTsize**2,1)
    end do
    endif
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
        call img2d_lowpass(a0,FFTsize,maxR)
        call img2d_norm(a0,FFTsize)
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if

        call img2d_mask(a0,FFTsize,imgmask,particle(i)%x,particle(i)%y)
        call img2d_norm(a0,FFTsize)
        ort60(1)=particle(i)%theta
        ort60(2)=particle(i)%phi
        ort60(3)=particle(i)%omega
        cent60(1)=-particle(i)%x
        cent60(2)=-particle(i)%y
        if(Ncycle .ne. 0) then
            do  nproc=1,Ncycle
                  dang =searchstep*(0.5)**(Nproc-1)
                  dtran=1.0*(0.5)**(Nproc-1)
                  call ort_cent_11(ort60(1),ort60(2),ort60(3),cent60(1),cent60(2),&
                  dang,dtran,ortcent)

                do jj=1,11
                     call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ortcent(1,jj),ortcent(2,jj)&
                     ,ortcent(3,jj),a1)
                     call img2d_mask(a1,FFTsize,imgmask,.0,.0)
                     call img2d_norm(a1,FFTsize)
                     call img2d_centshift_FT(a1,FFTsize,ortcent(4,jj),ortcent(5,jj))
                     call img2d_applyCTF(a1,FFTsize,ctf2d)
                     call img2d_lowpass(a1,FFTsize,maxR)
                     call img2d_norm(a1,FFTsize)


!            call cross_correlation(a2,a3,loc_box)
!                    Pr=a2(loc_box**2/2+loc_box/2+1)
                    Pr=phaseRedisual(a0,a1,FFTsize)
                    cent11(1,jj)=ortcent(4,jj)
                    cent11(2,jj)=ortcent(5,jj)
                    cent11(3,jj)=Pr
            end do
                   NNmax=maxloc(cent11(3,:))
                   ort60(1)=ortcent(1,NNmax(1))
                   ort60(2)=ortcent(2,NNmax(1))
                   ort60(3)=ortcent(3,NNmax(1))
                   cent60(1)=cent11(1,NNmax(1))
                   cent60(2)=cent11(2,NNmax(1))
                   cent60(3)=cent11(3,NNmax(1))
            end do
        end if
        particle(i)%dtheta=ort60( 1)-particle(i)%theta
        particle(i)%dphi=ort60( 2)-particle(i)%phi
        particle(i)%domega=ort60( 3)-particle(i)%omega
        particle(i)%dx=-cent60(1)-particle(i)%x
        particle(i)%dy=-cent60(2)-particle(i)%y

        particle(i)%theta=ort60( 1)
        particle(i)%phi  =ort60( 2)
        particle(i)%omega=ort60( 3)
        particle(i)%x    =-cent60(1)
        particle(i)%y    =-cent60(2)
        if(cent60(3).le.0) cent60(3)=.0
        if(cent60(3).ge.1) cent60(3)=1
        if(abs(cent60(1)).ge.FFTsize/2 .or.abs(cent60(2)).ge.FFTsize/2)then
            cent60(1)=0.0
            cent60(2)=0.0
            cent60(3)=0.0
        endif
        if(abs(ort60(1)).gt.999.0 .or.abs(ort60(2)).gt.999.0 .or.abs(ort60(3)).gt.999.0)then
            ort60(:)=0.0
            cent60(:)=0.0
        endif

        if(proc2d%pr_cc.eq.'y')then
           particle(i)%PR   =1.0-cent60(3)
        else
           particle(i)%PR= cent60(3)
        endif
        print*,i,cent60(:)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,particle(i)%astigang,particle(i)%dtheta,particle(i)%dphi,&
                      particle(i)%domega,particle(i)%dx,particle(i)%dy
    end do
    deallocate(a,a0,a1,a2,a3,loc_img2d,loc_proj2d,loc_img2d1,loc_proj2d1,loc_img2d2,loc_proj2d2)
    deallocate(core_FTr,core_FTi,a4)
    close(21)
    if(proc2d%imgstck.eq.'y') close(12)

900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
    return
end
