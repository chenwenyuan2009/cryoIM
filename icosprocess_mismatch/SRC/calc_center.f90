
subroutine calc_center
use common_icosproc
    implicit none
    integer i,j,ix,iy,iz,hx,ky,Nmax(1),ix0,iy0,imgID0
    real cc,z(-1:1,-1:1),dx,dy,centx,centy
    real,allocatable:: a1(:),a2(:),a(:)
    real,allocatable::img3d(:),img3di(:),shift3d(:)

    allocate(img3d(  FFTsize**3))
    allocate(img3di( FFTsize**3))
    allocate(shift3d(FFTsize**3))
    allocate(a1(FFTsize**2),a(FFTsize**2),a2(FFTsize**2))

!    open(21,file=trim(imgstck),form='unformatted',access='stream',status='old')
    open(31,file=trim(model3d),form='unformatted',access='stream',status='old')

    imgID0=-100
!    read(21) mrc
    read(31) mrc
    do i=1,FFTsize
      read(31) a
      do j=1,FFTsize**2
        img3d((i-1)*FFTsize**2+j)=a(j)
      end do
    end do
    close(31)

    do i=1,FFTsize**3
        iz=(i-1)/FFTsize**2
        iy=mod(i-1,FFTsize**2)/FFTsize
        ix=mod(mod(i-1,FFTsize**2),FFTsize)
        shift3d(i)=(-1)**(ix+iy+iz)
    end do

    if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) then
        call img3d_mask(img3d,FFTsize,imgmask,.0,.0,.0)
    end if
    img3d=img3d*shift3d
    img3di=.0
    call img3d_FFT(img3d,img3di,FFTsize,FFTsize,FFTsize,1)
    img3d =img3d *shift3d
    img3di=img3di*shift3d

    stckfile='abcdef'
    do i=first,last
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            close(11)
            stckfile=trim(particle(i)%stckfile)
            open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(11) mrc
            call fseek(11,4*FFTsize**2*(particle(i)%imgparticleID-1),1)
            read(11) a1
        else
            call fseek(11,4*FFTsize**2*(particle(i)%imgparticleID-particle(i-1)%imgparticleID-1),1)
            read(11) a1
        end if
        if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a1,FFTsize,imgmask,particle(i)%x,particle(i)%y)
        call img2d_lowpass(a1,FFTsize,maxR)

        call project2d_FT3dTo2d(img3d,img3di,FFTsize,particle(i)%theta,particle(i)%phi,particle(i)%omega,a2)
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
            call img2d_applyCTF(a2,FFTsize,ctf2d)
        end if
        if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a2,FFTsize,imgmask,.0,.0)
        call img2d_lowpass(a2,FFTsize,maxR)

        call img2d_norm(a1,FFTsize)
        call img2d_norm(a2,FFTsize)
        call cross_correlation(a1,a2,FFTsize)

        Nmax=maxloc(a1)-1
        iy0=Nmax(1)/FFTsize
        ix0=mod(Nmax(1),FFTsize)
        do ky=-1,1
            do hx=-1,1
                ix=ix0+hx
                iy=iy0+ky
                j=iy*FFTsize+ix+1
                z(hx,ky)=a1(j)
            end do
        end do
        call guass_regression(z,dx,dy,cc)
        centx=float(ix0-FFTsize/2)+dx
        centy=float(iy0-FFTsize/2)+dy
        print*,i,centx,centy,cc
        if(cc.ge. 1.0) cc=1.0
        if(cc.le.-1.0) cc=-1.0
        particle(i)%x0 =particle(i)%x
        particle(i)%y0 =particle(i)%y
        particle(i)%PR0=particle(i)%PR
        particle(i)%x  =centx
        particle(i)%y  =centy
        particle(i)%PR =acos(cc)*45.0/atan2(1.0,1.0)
        particle(i)%dx =particle(i)%x -particle(i)%x0
        particle(i)%dy =particle(i)%y -particle(i)%y0
        particle(i)%dPR=particle(i)%PR-particle(i)%PR0
    end do
    close(11)
    close(31)

    deallocate(img3d)
    deallocate(img3di)
    deallocate(shift3d)
    deallocate(a,a1,a2)
    return
 end subroutine






