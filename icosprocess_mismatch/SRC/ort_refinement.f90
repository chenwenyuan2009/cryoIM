subroutine ort_refinement
    use common_icosproc
    implicit none
    integer i,hx,ky,Nproc,jj,minN(1),imgID0
    real theta0,phi0,omega0,x0,y0,PR0
    real dang,dtran,ortcent(5,0:10)
    real PR(0:10),Fcut1(minR:maxR)
    real shift2d(FFTsize**2),a(FFTsize**2),b(FFTsize**2)
    real,allocatable::img2dstck(:)
    !open(31,file=trim(imgstck),form='unformatted',access='stream',status='old')

    imgID0=-100
    do i=1,FFTsize*FFTsize
        ky=(i-1)/FFTsize
        hx=mod(i-1,FFTsize)
        shift2d(i)=(-1)**(hx+ky)
    end do

    open(41,file=trim(newortfile))
    print*,trim(newortfile)
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))
    do i=first,last
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
!            close(11)
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(11) mrc
            allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
!            call fseek(11,4*FFTsize**2*(particle(i)%imgparticleID-1),1)
!            read(11) a
!        else
!            call fseek(11,4*FFTsize**2*(particle(i)%imgparticleID-particle(i-1)%imgparticleID-1),1)
            read(11) img2dstck
            close(11)
        end if
        a(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        a=a*shift2d
        b=.0
        if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a,FFTsize,imgmask,particle(i)%x,particle(i)%y)
        call img2d_FFT(a,b,FFTsize,FFTsize,1)
!        call img2d_FT_centshift(a,b,FFTsize,particle(i)%x0,particle(i)%y0)
        call FT_amp_expfit(a,b,FFTsize,minR,maxR,Fcut1)
        Fcut(:,0)=Fcut1(:)
        template2d(:,1,0)=a(:)
        template2d(:,2,0)=b(:)
        theta0=particle(i)%theta
        phi0  =particle(i)%phi
        omega0=particle(i)%omega
        x0    =particle(i)%x
        y0    =particle(i)%y

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if

        do Nproc=1,Ncycle
            dang =deltaAngle*(0.8)**(Nproc-1)
            dtran=deltaTrans*(0.8)**(Nproc-1)
            call ort_cent_11(theta0,phi0,omega0,x0,y0,dang,dtran,ortcent)
            do jj=0,10
                template(0)%theta=ortcent(1,jj)
                template(0)%phi  =ortcent(2,jj)
                template(0)%omega=ortcent(3,jj)
                template(0)%x    =ortcent(4,jj)
                template(0)%y    =ortcent(5,jj)
                call cal_crosslin_PR(PR(jj))
!                print*,i,Nproc,jj,dang,dtran
!                print*,ortcent(:,jj)
!                print*,PR(jj)
!                pause
            end do
            PR0 =minval(PR)
            minN=minloc(PR)
            theta0=ortcent(1,minN(1)-1)
            phi0  =ortcent(2,minN(1)-1)
            omega0=ortcent(3,minN(1)-1)
            x0    =ortcent(4,minN(1)-1)
            y0    =ortcent(5,minN(1)-1)
!            print*,i,Nproc,minN,PR0
!            print*,theta0,phi0,omega0,x0,y0
!            pause
        end do
        particle(i)%theta0=particle(i)%theta
        particle(i)%phi0  =particle(i)%phi
        particle(i)%omega0=particle(i)%omega
        particle(i)%x0    =particle(i)%x
        particle(i)%y0    =particle(i)%y
        particle(i)%PR0   =particle(i)%PR

        particle(i)%theta=theta0
        particle(i)%phi  =phi0
        particle(i)%omega=omega0
        particle(i)%x    =x0
        particle(i)%y    =y0
        particle(i)%PR   =PR0

        particle(i)%dtheta=particle(i)%theta-particle(i)%theta0
        particle(i)%dphi  =particle(i)%phi  -particle(i)%phi0
        particle(i)%domega=particle(i)%omega-particle(i)%omega0
        particle(i)%dx    =particle(i)%x    -particle(i)%x0
        particle(i)%dy    =particle(i)%y    -particle(i)%y0
        particle(i)%dPR   =particle(i)%PR   -particle(i)%PR0
        print*,i,theta0,phi0,omega0,x0,y0,PR0

        write(41,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
        particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1,       particle(i)%df2,   particle(i)%astigang, &
                      particle(i)%dtheta,particle(i)%dphi,particle(i)%domega, &
                      particle(i)%dx,particle(i)%dy,particle(i)%dPR

    end do
    close(11)
    close(41)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
    return
end
