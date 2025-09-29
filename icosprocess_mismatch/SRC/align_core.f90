subroutine algin_core
use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1,NNmax(1)
    integer ix,iy,n1,nproc,jj,ix0,iy0
    real dx,dy,cc,dang,dtran,ortcent(5,7)
    real ort60(3,60),cent60(3,60),Cmax,z(-1:1,-1:1),cent11(3,7)
    real(kind=8) Pr,phaseRedisual
    real,allocatable::a0(:),a1(:),a2(:),a3(:),img2dstck(:)
    real,allocatable::core_FTr(:),core_FTi(:)


    allocate(a0(FFTsize**2),a1(FFTsize**2),a2(FFTsize**2),img2dstck(FFTsize**2))
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),a3(FFTsize**2))

    print*,'reading 3D map...'
    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    if(proc2d%imgstck.eq.'y')then
            open(12,file=trim(imgstck),form='unformatted',access='stream',status='old')
            read(12) mrc
    endif
    open(21,file=trim(newortfile))
    read(11) mrc
    do i=1,FFTsize
      read(11) a3
      do j=1,FFTsize**2
        core_FTr((i-1)*FFTsize**2+j)=a3(j)
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
        call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)
!        particle(i)%x=.0
!        particle(i)%y=.0
        call img2d_lowpass(a0,FFTsize,maxR)
        call img2d_norm(a0,FFTsize)
        call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        call img2d_norm(a0,FFTsize)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if

        call icos_equalort(particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60)
        !print*,particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60(:,1)
        do j=1,60
               call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ort60(1,j),ort60(2,j),ort60(3,j),a1)
               call img2d_mask(a1,FFTsize,imgmask,.0,.0)
               call img2d_norm(a1,FFTsize)
               call img2d_applyCTF(a1,FFTsize,ctf2d)
               call img2d_lowpass(a1,FFTsize,maxR)
               call img2d_norm(a1,FFTsize)
               a2=a0
               !call cross_correlation(a2,a1,FFTsize)
               !Pr=a2(FFTsize**2/2+FFTsize/2+1)
               Pr=phaseRedisual(a2,a1,FFTsize)

               if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) then
                Nmax=maxloc(a2)-1
                iy0=Nmax(1)/FFTsize
                ix0=mod(Nmax(1),FFTsize)
                 if((iy0 .lt. (FFTsize-2)).and.(ix0 .lt.(FFTsize-2)))then
                    if((iy0.gt.2).and.(ix0.gt.2))then
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
                    cent60(1,j)=.0
                    cent60(2,j)=.0
                    cent60(3,j)=.0
                    endif
                    else
                    cent60(1,j)=.0
                    cent60(2,j)=.0
                    cent60(3,j)=.0
                endif
!                cc=maxval(a2)
!                cent60(3,j)=cc
!                cent60(1,j)=float(ix0-FFTsize/2)
!                cent60(2,j)=float(iy0-FFTsize/2)
               else
                cent60(1,j)=particle(i)%x
                cent60(2,j)=particle(i)%y
                cent60(3,j)=Pr!a2(FFTsize**2/2+FFTsize/2+1)
               end if
 !           end if
        end do
        Nmax=maxloc(cent60(3,:))
        if(cent60(3,Nmax(1)).le.0) cent60(3,Nmax(1))=.0
        if(cent60(3,Nmax(1)).ge.1) cent60(3,Nmax(1))=1
        if(abs(cent60(1,Nmax(1))).ge.FFTsize/2 .or.abs(cent60(2,Nmax(1))).ge.FFTsize/2)then
            cent60(1,Nmax(1))=0.0
            cent60(2,Nmax(1))=0.0
            cent60(3,Nmax(1))=0.0
        endif
        particle(i)%theta=ort60( 1,Nmax(1))
        particle(i)%phi  =ort60( 2,Nmax(1))
        particle(i)%omega=ort60( 3,Nmax(1))
        particle(i)%x    =cent60(1,Nmax(1))
        particle(i)%y    =cent60(2,Nmax(1))
        if(proc2d%pr_cc.eq.'y')then
           particle(i)%PR   =1.0-cent60(3,Nmax(1))
        else
           particle(i)%PR= cent60(3,Nmax(1))
        endif

        print*,i,Nmax,cent60(:,Nmax(1))

        Cmax=acos(Cmax)*45.0/atan2(1.0,1.0)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,   particle(i)%astigang,Nmax!,cent60(3,1),cent60(3,2),cent60(3,3),&
                      !cent60(3,4),cent60(3,5)
    end do
    deallocate(a3,a0,a1,a2)
    deallocate(core_FTr,core_FTi)
    close(21)
    if(proc2d%imgstck.eq.'y') close(12)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5)!,f10.4,f10.4,f10.4,f10.4,f10.4
    return
end
