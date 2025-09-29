
subroutine ort_globesearch
    use common_icosproc
    implicit none
    integer i,hx,ky,minN(3)
    integer itheta,iphi,Ntheta,Nphi,Nomega,imgID0
    real theta,phi,omega,PR0,tmp
    real deg2rad,A1,B1,C1,A2,B2,C2,x,y,z,z0,z1,z2

    real,allocatable:: Fcut1(:),shift2d(:),a(:),b(:)
    real,allocatable:: PR_array(:,:,:),PR_omega(:),img2dstck(:)

    allocate(Fcut1(minR:maxR),shift2d(FFTsize**2),a(FFTsize**2),b(FFTsize**2))

    imgID0=-100
    deg2rad=atan2(1.0,1.0)/45.0
    A1=-sin(31.717*deg2rad)*cos(69.09*deg2rad);B1=-cos(31.717*deg2rad)*cos(69.09*deg2rad);C1=sin(31.717*deg2rad)*sin(69.09*deg2rad)
    A2=-sin(31.717*deg2rad)*cos(69.09*deg2rad);B2= cos(31.717*deg2rad)*cos(69.09*deg2rad);C2=sin(31.717*deg2rad)*sin(69.09*deg2rad)
    Ntheta=int(21.0 /searchstep)
    Nphi  =int(32.0 /searchstep)
    Nomega=int(360.0/searchstep)
    allocate(PR_array(0:Nomega-1,-Nphi:Nphi,0:Ntheta))
    allocate(PR_omega(0:Nomega-1))
    PR_array=1000.0
    do i=1,FFTsize*FFTsize
        ky=(i-1)/FFTsize
        hx=mod(i-1,FFTsize)
        shift2d(i)=(-1)**(hx+ky)
    end do

    stckfile='abcdef'
    open(41,file=trim(newortfile))
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
        call img2d_mask(a,FFTsize,imgmask,particle(i)%x,particle(i)%y)
        a=a*shift2d
        b=.0
        call img2d_FFT(a,b,FFTsize,FFTsize,1)
        call FT_amp_expfit(a,b,FFTsize,minR,maxR,Fcut1)  !fit the structure factor
        call img2d_FT_centshift(a,b,FFTsize,particle(i)%x,particle(i)%y)  !shift center in Fourier space
        Fcut(:,0)=Fcut1(:)
        template2d(:,1,0)=a(:)
        template2d(:,2,0)=b(:)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if

        do itheta=0,Ntheta
!            print*,'itheta=',itheta
            theta=69.0+searchstep*itheta
            tmp=theta*deg2rad
            z= cos(tmp)
            z0=sin(tmp)
            do iphi=0,Nphi
                phi=searchstep*iphi
!                print*,'iphi=',iphi,theta,phi
                tmp=phi*deg2rad
                x=z0*cos(tmp)
                y=z0*sin(tmp)
                z1=-(A1*x+B1*y)/C1
                z2=-(A2*x+B2*y)/C2
                if((z.ge..0).and.(z.le.z1).and.(z.le.z2)) then
                    omega=.0
!                    print*,theta,phi,omega
                    template(0)%theta=theta
                    template(0)%phi  =phi
                    template(0)%omega=.0
                    template(0)%x    =particle(i)%x
                    template(0)%y    =particle(i)%y
                    call cal_crosslin_PR_omega(Nomega,PR_omega)
                    PR_array(:,iphi,itheta)=PR_omega(:)
                    template(0)%theta=theta
                    template(0)%phi  =-phi
                    template(0)%omega=.0
                    template(0)%x    =particle(i)%x
                    template(0)%y    =particle(i)%y
                    call cal_crosslin_PR_omega(Nomega,PR_omega)
                    PR_array(:,-iphi,itheta)=PR_omega(:)
                end if
            end do
        end do
        PR0 =minval(PR_array)
        minN=minloc(PR_array)
        theta=69.0+searchstep*(minN(3)-1)
        phi  =searchstep*(minN(2)-Nphi-1)
        omega=searchstep*(minN(1)-1)
        print*,i,theta,phi,omega,particle(i)%x,particle(i)%y,PR0

        particle(i)%theta0=particle(i)%theta
        particle(i)%phi0  =particle(i)%phi
        particle(i)%omega0=particle(i)%omega
        particle(i)%x0    =particle(i)%x
        particle(i)%y0    =particle(i)%y
        particle(i)%PR0   =particle(i)%PR

        particle(i)%theta=theta
        particle(i)%phi  =phi
        particle(i)%omega=omega
        particle(i)%PR   =PR0

        particle(i)%dtheta=particle(i)%theta-particle(i)%theta0
        particle(i)%dphi  =particle(i)%phi  -particle(i)%phi0
        particle(i)%domega=particle(i)%omega-particle(i)%omega0
        particle(i)%dx    =particle(i)%x    -particle(i)%x0
        particle(i)%dy    =particle(i)%y    -particle(i)%y0
        particle(i)%dPR   =particle(i)%PR   -particle(i)%PR0

        write(41,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
        particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1,       particle(i)%df2,   particle(i)%astigang, &
                      particle(i)%dtheta,particle(i)%dphi,particle(i)%domega, &
                      particle(i)%dx,particle(i)%dy,particle(i)%dPR

    end do
    deallocate(PR_array,Fcut1,shift2d,a,b)
    deallocate(PR_omega)
    close(11)
    close(41)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
    return
end

