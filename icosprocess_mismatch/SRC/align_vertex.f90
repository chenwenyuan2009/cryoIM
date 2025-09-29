subroutine algin_vertex
use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1,NNmax(1),nproc
    integer ix0,iy0,ix,iy,n1,ip1,ip2,x0,y0,ix00,iy00,jj,ii
    real dx,dy,cc,dn,dm,RotMat(3,3),delts,dang,dtran,ortcent(5,11),dfcent(3,7),dfpr(7)
    real ort60(3,5),cent60(3,5),Cmax,z(-1:1,-1:1),initiald(3,1),newd(3,1),cent11(3,11)
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
    !open(20,file=trim(result3d),form='unformatted',access='stream')
    open(21,file=trim(newortfile))
    read(11) mrc
    do i=1,FFTsize
      read(11) a
      do j=1,FFTsize**2
        core_FTr((i-1)*FFTsize**2+j)=a(j)
      end do
    end do

    close(11)
!    mrc%nx=loc_box
!    mrc%ny=loc_box
!    write(20) mrc
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
    stckfile='abcd'
    !read(12) mrc
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
        if(proc2d%centshift.eq.'y') then
        call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)
        particle(i)%x=.0
        particle(i)%y=.0
        endif
        call img2d_lowpass(a0,FFTsize,maxR)
        !call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        call img2d_norm(a0,FFTsize)
        a4=a0
        call img2d_centshift_FT(a4,FFTsize,particle(i)%x,particle(i)%y)
        call img2d_lowpass(a4,FFTsize,maxR)
        !call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        call img2d_norm(a4,FFTsize)
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
        call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if
        call icos_equalort05(particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60)
        !print*,particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60(:,1)
        do j=1,5
            call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ort60(1,j),ort60(2,j),ort60(3,j),a1)
            call img2d_applyCTF(a1,FFTsize,ctf2d)
            call img2d_lowpass(a1,FFTsize,maxR)
            !call img2d_mask(a1,FFTsize,imgmask,.0,.0)
            call img2d_norm(a1,FFTsize)

            !ort60(2,j)=ort60(2,j)
            !ort60(3,j)=ort60(3,j)
            call euler_mat_vertex(ort60(1,j),ort60(2,j),ort60(3,j),RotMat)
             newd=matmul(RotMat,initiald)
            x0=int(newd(1,1))+FFTsize/2
            y0=int(newd(2,1))+FFTsize/2
            ix00=int(newd(1,1))
            iy00=int(newd(2,1))

            do iy=y0-loc_box/2,y0+loc_box/2-1
                do ix=x0-loc_box/2,x0+loc_box/2-1
                    ip1=iy*FFTsize+ix+1
                    ip2=(iy-y0+loc_box/2)*loc_box+ix-x0+loc_box/2+1
                    loc_img2d(ip2,j) =a4(ip1)
                    loc_proj2d(ip2,j)=a1(ip1)
!                    if(ix00 .ne. 0)then
!                    delts=abs((iy00/ix00)*(ix-ix00)-(iy-iy00))/sqrt(float((iy00/ix00))**2-1.0)
!                    if(delts.gt.50.0)then
!                        loc_img2d(ip2,j)=.0
!                        loc_proj2d(ip2,j)=.0
!                    endif
!                    endif
                end do
            end do

            a2=loc_img2d(:,j)
            a3=loc_proj2d(:,j)
            call img2d_mask(a2,loc_box,imgmask,.0,.0)
            call img2d_mask(a3,loc_box,imgmask,.0,.0)
            call img2d_norm(a2,loc_box)
            call img2d_norm(a3,loc_box)

!            call cross_correlation(a2,a3,loc_box)
!            Pr=a2(loc_box**2/2+loc_box/2+1)
            Pr=phaseRedisual(a2,a3,loc_box)
            cent60(3,j)=Pr
        end do
        Nmax=maxloc(cent60(3,:))
        cent60(1,Nmax(1))=-particle(i)%x
        cent60(2,Nmax(1))=-particle(i)%y
        if(Ncycle .ne. 0) then
           cent60(1,Nmax(1))=-particle(i)%x
           cent60(2,Nmax(1))=-particle(i)%y
            do  nproc=1,Ncycle
                  dang =searchstep*(0.8)**(Nproc-1)
                  dtran=2.0*(0.8)**(Nproc-1)
                  call ort_cent_11(ort60(1,Nmax(1)),ort60(2,Nmax(1)),ort60(3,Nmax(1)),cent60(1,Nmax(1)),cent60(2,Nmax(1)),&
                  dang,dtran,ortcent)

                do jj=1,11
                    !print*,ortcent(:,jj)
                     call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ortcent(1,jj),ortcent(2,jj)&
                     ,ortcent(3,jj),a1)
                     call img2d_centshift_FT(a1,FFTsize,ortcent(4,jj),ortcent(5,jj))
                     call img2d_applyCTF(a1,FFTsize,ctf2d)
                     call img2d_lowpass(a1,FFTsize,maxR)
                     !call img2d_mask(a1,FFTsize,imgmask,.0,.0)
                     call img2d_norm(a1,FFTsize)
            call euler_mat_vertex(ortcent(1,jj),ortcent(2,jj),ortcent(3,jj),RotMat)
            newd=matmul(RotMat,initiald)

            x0=int(newd(1,1))+FFTsize/2
            y0=int(newd(2,1))+FFTsize/2
            ix00=int(newd(1,1))
            iy00=int(newd(2,1))

            do iy=y0-loc_box/2,y0+loc_box/2-1
                do ix=x0-loc_box/2,x0+loc_box/2-1
                    ip1=iy*FFTsize+ix+1
                    ip2=(iy-y0+loc_box/2)*loc_box+ix-x0+loc_box/2+1
                    loc_img2d1(ip2,jj) =a0(ip1)
                    loc_proj2d1(ip2,jj)=a1(ip1)
!                    if(ix00 .ne. 0)then
!                    delts=abs((iy00/ix00)*(ix-ix0)-(iy-iy0))/sqrt(float((iy00/ix00))**2-1.0)
!                    if(delts.lt.50.0)then
!                        loc_img2d(ip2,j)=.0
!                        loc_proj2d(ip2,j)=.0
!                    endif
!                    endif
                end do
            end do

            a2=loc_img2d1(:,jj)
            a3=loc_proj2d1(:,jj)
            call img2d_mask(a2,loc_box,imgmask,.0,.0)
            call img2d_mask(a3,loc_box,imgmask,.0,.0)
            call img2d_norm(a2,loc_box)
            call img2d_norm(a3,loc_box)

!            call cross_correlation(a2,a3,loc_box)
!                    Pr=a2(loc_box**2/2+loc_box/2+1)
                    Pr=phaseRedisual(a2,a3,loc_box)
                    cent11(1,jj)=ortcent(4,jj)
                    cent11(2,jj)=ortcent(5,jj)
                    cent11(3,jj)=Pr
            end do
                   NNmax=maxloc(cent11(3,:))
                   ort60(1,Nmax(1))=ortcent(1,NNmax(1))
                   ort60(2,Nmax(1))=ortcent(2,NNmax(1))
                   ort60(3,Nmax(1))=ortcent(3,NNmax(1))
                   cent60(1,Nmax(1))=cent11(1,NNmax(1))
                   cent60(2,Nmax(1))=cent11(2,NNmax(1))
                   cent60(3,Nmax(1))=cent11(3,NNmax(1))
            end do

        end if
        particle(i)%dtheta=ort60( 1,Nmax(1))-particle(i)%theta
        particle(i)%dphi=ort60( 2,Nmax(1))-particle(i)%phi
        particle(i)%domega=ort60( 3,Nmax(1))-particle(i)%omega
        particle(i)%dx=-cent60(1,Nmax(1))-particle(i)%x
        particle(i)%dy=-cent60(2,Nmax(1))-particle(i)%y

        particle(i)%theta=ort60( 1,Nmax(1))
        particle(i)%phi  =ort60( 2,Nmax(1))
        particle(i)%omega=ort60( 3,Nmax(1))
        particle(i)%x    =-cent60(1,Nmax(1))
        particle(i)%y    =-cent60(2,Nmax(1))
!        if(Ncycle .ne. 0) then
!          dforigin(1)=particle(i)%df1
!          dforigin(2)=particle(i)%df2
!          dforigin(3)=particle(i)%astigang
!          call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)
!          call img2d_lowpass(a0,FFTsize,maxR)
!          call img2d_norm(a0,FFTsize)
!          call img2d_mask(a0,FFTsize,imgmask,.0,.0)
!          call img2d_norm(a0,FFTsize)
!          call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,particle(i)%theta,particle(i)%phi&
!               ,particle(i)%omega,a1)
!          call img2d_lowpass(a1,FFTsize,maxR)
!          call img2d_mask(a1,FFTsize,imgmask,.0,.0)
!          call img2d_norm(a1,FFTsize)
!          do  nproc=1,Ncycle
!              dang =10*(0.8)**(Nproc-1)
!              dtran=0.05*(0.8)**(Nproc-1)
!          call dfrefine(dforigin(1),dforigin(2),dforigin(3),dtran,dang,dfcent)
!          do ii=1,7
!           call getCTF2d(FFTsize,dfcent(1,ii),dfcent(2,ii),dfcent(3,ii),apix,&
!                CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
!           call img2d_applyCTF(a1,FFTsize,ctf2d)
!           call img2d_norm(a1,FFTsize)
!
!!            call cross_correlation(a2,a3,loc_box)
!!            Pr=a2(loc_box**2/2+loc_box/2+1)
!            Pr=phaseRedisual(a0,a1,FFTsize)
!            dfpr(ii)=Pr
!          end do
!          NNmax=maxloc(dfpr(:))
!          dforigin(1)=dfcent(1,NNmax(1))
!          dforigin(2)=dfcent(2,NNmax(1))
!          dforigin(3)=dfcent(3,NNmax(1))
!          cent60(3,Nmax(1))=dfpr(NNmax(1))
!          end do
!          particle(i)%deltdf1=dforigin(1)-particle(i)%df1
!          particle(i)%deltdf2=dforigin(2)-particle(i)%df2
!          particle(i)%deltdfa=dforigin(3)-particle(i)%astigang
!          particle(i)%df1=dforigin(1)
!          particle(i)%df2=dforigin(2)
!          particle(i)%astigang=dforigin(3)
!        endif
        if(cent60(3,Nmax(1)).le.0) cent60(3,Nmax(1))=.0
        if(cent60(3,Nmax(1)).ge.1) cent60(3,Nmax(1))=1
        if(abs(cent60(1,Nmax(1))).ge.FFTsize/2 .or.abs(cent60(2,Nmax(1))).ge.FFTsize/2)then
            cent60(1,Nmax(1))=0.0
            cent60(2,Nmax(1))=0.0
            cent60(3,Nmax(1))=0.0
        endif
        if(abs(ort60(1,Nmax(1))).gt.999.0 .or.abs(ort60(2,Nmax(1))).gt.999.0 .or.abs(ort60(3,Nmax(1))).gt.999.0)then
            ort60(:,Nmax(1))=0.0
            cent60(:,Nmax(1))=0.0
        endif

        if(proc2d%pr_cc.eq.'y')then
           particle(i)%PR   =1.0-cent60(3,Nmax(1))
        else
           particle(i)%PR= cent60(3,Nmax(1))
        endif
        print*,i,Nmax,cent60(:,Nmax(1))
        !write(20)loc_proj2d(:,Nmax(1))

 !       Cmax=acos(Cmax)*45.0/atan2(1.0,1.0)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,particle(i)%astigang,Nmax,particle(i)%dtheta,particle(i)%dphi,&
                      particle(i)%domega,particle(i)%dx,particle(i)%dy!,particle(i)%deltdf1,particle(i)%deltdf2,&
                      !particle(i)%deltdfa

    end do
    deallocate(a,a0,a1,a2,a3,loc_img2d,loc_proj2d,loc_img2d1,loc_proj2d1,loc_img2d2,loc_proj2d2)
    deallocate(core_FTr,core_FTi,a4)
    close(21)
    if(proc2d%imgstck.eq.'y') close(12)
    !close(20)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5,f10.4,f10.4,f10.4,f10.4,f10.4)!,&
          ! f10.4,f10.4,f10.4)
    return
end
