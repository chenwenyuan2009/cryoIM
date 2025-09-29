subroutine loc_search
use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1,NNmax(1)
    integer ix0,iy0,ix,iy,n1,nproc,jj
    real dx,dy,cc,dang,dtran,ortcent(5,7)
    real ort60(3,60),cent60(3,60),Cmax,z(-1:1,-1:1),cent11(3,7)
    real,allocatable::a0(:),a1(:),a2(:),a3(:)
    real,allocatable::core_FTr(:),core_FTi(:)


    allocate(a0(FFTsize**2),a1(FFTsize**2),a2(FFTsize**2))
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),a3(FFTsize**2))

    print*,'reading 3D map...'
    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    open(12,file=trim(imgstck),form='unformatted',access='stream',status='old')
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
    read(12) mrc
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
        read(12) a0
        call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)
        particle(i)%x=.0
        particle(i)%y=.0
        call img2d_lowpass(a0,FFTsize,maxR)
        call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        call img2d_norm(a0,FFTsize)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if
        if(Ncycle .ne. 0) then
            do  nproc=1,Ncycle
                  dang =searchstep*(0.8)**(Nproc-1)
                  dtran=searchstep*(0.8)**(Nproc-1)
                  call ort_cent_11(particle(i)%theta,particle(i)%phi,particle(i)%omega,particle(i)%x,particle(i)%y,&
                  dang,dtran,ortcent)
                  do jj=0,6
                     call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ortcent(1,jj),ortcent(2,jj)&
                     ,ortcent(3,jj),a1)
                     call img2d_applyCTF(a1,FFTsize,ctf2d)
                     call img2d_lowpass(a1,FFTsize,maxR)
                     call img2d_mask(a1,FFTsize,imgmask,.0,.0)
                     call img2d_norm(a1,FFTsize)
                     a2=a0
                     call cross_correlation(a2,a1,FFTsize)
                     if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) then
                        NNmax=maxloc(a2)-1
                        iy0=NNmax(1)/FFTsize
                        ix0=mod(NNmax(1),FFTsize)
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
                          cent11(1,jj)=float(ix0-FFTsize/2)+dx
                          cent11(2,jj)=float(iy0-FFTsize/2)+dy
                          cent11(3,jj)=cc
                          else
                          cent11(1,jj)=.0
                          cent11(2,jj)=.0
                          cent11(3,jj)=a2(FFTsize**2/2+FFTsize/2+1)
                          endif
                         else
                          cent11(1,jj)=.0
                          cent11(2,jj)=.0
                          cent11(3,jj)=a2(FFTsize**2/2+FFTsize/2+1)
                         endif
!                          cc=maxval(a2)
!                          cent11(3,jj)=cc
!                          cent11(1,jj)=float(ix0-FFTsize/2)
!                          cent11(2,jj)=float(iy0-FFTsize/2)
                      else
                        cent11(1,jj)=.0
                        cent11(2,jj)=.0
                        cent11(3,jj)=a2(FFTsize**2/2+FFTsize/2+1)
                      end if
                   end do
                   NNmax=maxloc(cent11(3,:))
                   particle(i)%theta=ortcent(1,NNmax(1))
                   particle(i)%phi  =ortcent(2,NNmax(1))
                   particle(i)%omega=ortcent(3,NNmax(1))
                   particle(i)%x    =cent11(1,NNmax(1))
                   particle(i)%y    =cent11(2,NNmax(1))
                   particle(i)%PR   =cent11(3,NNmax(1))
            end do
        end if

        print*,i,particle(i)%x,particle(i)%y,particle(i)%PR

        Cmax=acos(Cmax)*45.0/atan2(1.0,1.0)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,   particle(i)%astigang

    end do
    deallocate(a3,a0,a1,a2)
    deallocate(core_FTr,core_FTi)
    close(21)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5)
    return
end
