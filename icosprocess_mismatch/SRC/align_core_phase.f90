

subroutine algin_core_phase
use common_icosproc
    implicit none
    integer pID,Nmax(1),item,hx,ky,lz,i2,i3,n1,n2,imgID0,no0
    real theta,phi,omega,ort60(3,5),cc60(5),cent(2,5)
    real,allocatable:: a1(:),a2(:),a20(:)
    real,allocatable:: b1(:),b2(:)
    real cc,centx,centy,dr(4),x,y,xx,lowpass,Huv
    real,allocatable:: ft3dr(:),ft3di(:),img2dstck(:)
    real,allocatable:: shift2d(:)

    if(proc2d%imgstck.eq.'y')then
            open(31,file=trim(imgstck),form='unformatted',access='stream',status='old')
            read(31) mrc
    endif
    open(41,file=trim(newortfile))
    allocate(a1(FFTsize**2),a2(FFTsize**2),a20(FFTsize**2),img2dstck(FFTsize**2))
    allocate(b1(FFTsize**2),b2(FFTsize**2),shift2d(FFTsize**2))


    n1=FFTsize
    n2=FFTsize
    lowpass=2.0*maxR**2
    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        shift2d(i2)=(-1)**(hx+ky)
    end do

!********************************************
    open(21,file=trim(model3d),form='unformatted',access='stream',status='old')
    allocate(ft3di(FFTsize**3),ft3dr(FFTsize**3))
    read(21) mrc
    read(21) ft3dr
    close(21)
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            ft3dr(i3)=-ft3dr(i3)
        end if
    end do
    ft3di=.0
    stckfile='abcd'
    print*,'proferming FFT3d for 3d model'
    call img3d_FFT(ft3dr,ft3di,FFTsize,FFTsize,FFTsize,1)
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            ft3dr(i3)=-ft3dr(i3)
            ft3di(i3)=-ft3di(i3)
        end if
    end do

    !read(31) mrc
    if(proc2d%imgstck.eq.'y')then
    do pID=1,first-1
        call fseek(31,4*FFTsize**2,1)
    end do
    endif
    no0=-1
    imgID0=-100
    do pID=first,last
        theta=particle(pID)%theta
        phi  =particle(pID)%phi
        omega=particle(pID)%omega
    if(proc2d%imgstck.eq.'y')then
          read(31)a20
       else
        if(trim(stckfile).ne.trim(particle(pID)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(pID)%stckfile)
            open(31,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(31) mrc
            allocate(img2dstck(FFTsize**2*mrc%nz))
            read(31)img2dstck
            close(31)
        end if
        a20(:)=img2dstck((particle(pID)%imgparticleID-1)*FFTsize**2+1:particle(pID)%imgparticleID*FFTsize**2)
        endif
        call img2d_mask(a20,FFTsize,imgmask,particle(pID)%x,particle(pID)%y)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(pID)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(pID)%df1,particle(pID)%df2,particle(pID)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(pID)%imgID
            end if
        end if

        call icos_equalort05(theta,phi,omega,ort60)
        do item=1,5
            call project2d_FT3dTo2d(ft3dr,ft3di,FFTsize,ort60(1,item),ort60(2,item),ort60(3,item),a1)
            call img2d_applyCTF(a1,FFTsize,ctf2d)
            call img2d_mask(a1,FFTsize,imgmask,.0,.0)
            b1=.0
            a2=a20
            b2=.0
            call FFT1d(a1,b1,n1*n2,n1,n1,1)
            call FFT1d(a1,b1,n1*n2,n2,n1*n2,1)
            call FFT1d(a2,b2,n1*n2,n1,n1,1)
            call FFT1d(a2,b2,n1*n2,n2,n1*n2,1)
            xx=.0

            do i2=1,FFTsize**2
                ky=(i2-1)/FFTsize
                hx=mod(i2-1,FFTsize)
                dr(1)=float(hx**2+ky**2)
                dr(2)=float((hx-FFTsize)**2+(ky**2))
                dr(3)=float(hx**2+(ky-FFTsize)**2)
                dr(4)=float((hx-FFTsize)**2+(ky-FFTsize)**2)
                Huv=exp(-minval(dr)/lowpass)
                x=atan2(b1(i2),a1(i2))
                y=atan2(b2(i2),a2(i2))
                a1(i2)=cos(y-x)*Huv*(-1)**(hx+ky)
                b1(i2)=sin(y-x)*Huv*(-1)**(hx+ky)
                xx=xx+Huv
            end do

            call FFT1d(a1,b1,n1*n2,n1,n1,-1)
            call FFT1d(a1,b1,n1*n2,n2,n1*n2,-1)
            a1=a1/xx
            if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) then
                Nmax=maxloc(a1)-1
                centy=Nmax(1)/FFTsize-FFTsize/2
                centx=mod(Nmax(1),FFTsize)-FFTsize/2
                cc=maxval(a1)
            else
                centy=0
                centx=0
                cc=a1(FFTsize**2/2+FFTsize/2+1)
            end if
            CC60(item)=cc
            cent(1,item)=centx
            cent(2,item)=centy
        end do
        Nmax=maxloc(cc60)
        particle(pID)%theta=ort60(1,Nmax(1))
        particle(pID)%phi  =ort60(2,Nmax(1))
        particle(pID)%omega=ort60(3,Nmax(1))
        particle(pID)%x=cent(1,Nmax(1))
        particle(pID)%y=cent(2,Nmax(1))
        if(proc2d%pr_cc.eq.'y')then
           particle(pID)%PR   =1.0-maxval(cc60)
        else
           particle(pID)%PR= maxval(cc60)
        endif
        print*,pID,Nmax,cent(1,Nmax(1)),cent(2,Nmax(1)),particle(pID)%PR  !,ccavg,ccdev

        write(41,900) particle(pID)%stckfile,particle(pID)%particleID,particle(pID)%theta,particle(pID)%phi,particle(pID)%omega, &
                      particle(pID)%x,particle(pID)%y,particle(pID)%PR,particle(pID)%imgID,particle(pID)%imgparticleID,  &
                      particle(pID)%df1,       particle(pID)%df2,   particle(pID)%astigang,Nmax!,CC60(1),CC60(2),CC60(3),&
                     ! CC60(4),CC60(5)
    end do

    deallocate(a1,a2,a20,b1,b2)
    deallocate(ft3di,ft3dr)
    close(21)
    if(proc2d%imgstck.eq.'y')close(31)
    close(41)
    close(51)

900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5)!,f10.4,f10.4,f10.4,f10.4,f10.4
    return
 end subroutine








