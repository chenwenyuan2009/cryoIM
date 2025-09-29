subroutine phase_redsiual
    use common_icosproc
    implicit none
    integer i,hx,hy,jj,maxN(1),FR,i3,imgID0,j,ky,lz,no0,no1
    real theta0,phi0,omega0,x0,y0,z0,pi,Tpi
    real ort60(3,5),mincc,PR1,PR2,phas1,phas2,amp1,amp2,diff
    real pr(5),Fcut1(minR:maxR),Fcut2(minR:maxR),alpha1,alpha2
    real FT_commlin1(2,minR:maxR),FT_commlin2(2,minR:maxR)
    real,allocatable:: img2d(:),shift2d(:),a(:),b(:),a1(:),b1(:),a3(:)
    real,allocatable::core_FTr(:),core_FTi(:)
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),img2d(FFTsize**2),b1(FFTsize**2))
    allocate(shift2d(FFTsize**2),a(FFTsize**2),b(FFTsize**2),a1(FFTsize**2),a3(FFTsize**2))

    pi=4.0*atan2(1.0,1.0)
    Tpi=2.0*pi
    PR1=.0
    PR2=.0
    pr=.0
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
    do i=1,FFTsize**2
        ky=(i-1)/FFTsize
        hx=mod(i-1,FFTsize)
        shift2d(i)=(-1)**(hx+ky)
    end do

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
        read(12) img2d
        a=img2d*shift2d
        b=.0
        if(proc2d%centshift.eq.'y')then
            call img2d_centshift_FT(a,FFTsize,particle(i)%x,particle(i)%y)
            particle(i)%x=.0
            particle(i)%y=.0
        endif
        if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y'))call img2d_mask(a,FFTsize,imgmask,particle(i)%x,particle(i)%y)
        call img2d_FFT(a,b,FFTsize,FFTsize,1)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if

        call icos_equalort05(particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60)
        do j=1,5
               call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,ort60(1,j),ort60(2,j),ort60(3,j),a1)
               b1=.0
               call img2d_mask(a1,FFTsize,imgmask,.0,.0)
               call img2d_applyCTF(a1,FFTsize,ctf2d)
               call img2d_FFT(a1,b1,FFTsize,FFTsize,1)
               call crosslin_angle(particle(i)%theta,particle(i)%phi,particle(i)%omega,ort60(1,j),&
                    ort60(2,j),ort60(3,j),alpha1,alpha2)

               call img2d_fill(a,b,FFTsize,minR,maxR,alpha1,FT_commlin1)
               call img2d_fill(a1,b1,FFTsize,minR,maxR,alpha2,FT_commlin2)
            if((abs(alpha1).lt.999.0).and.(abs(alpha2).lt.999.0))then
               do FR=minR,maxR

                    amp1=sqrt(FT_commlin1(1,FR)**2+FT_commlin1(2,FR)**2)
                    amp2=sqrt(FT_commlin2(1,FR)**2+FT_commlin2(2,FR)**2)
                    if(diff.ge.pi) diff=Tpi-diff
                    phas1=atan2(FT_commlin1(2,FR),FT_commlin1(1,FR))
                    phas2=atan2(FT_commlin2(2,FR),FT_commlin2(1,FR))
                    diff=abs(phas2-phas1)
                    diff=amod(diff,Tpi)
                    amp1=0.5*(amp1+amp2)
                    PR1=PR1+diff*amp1
                    PR2=PR2+amp1
                end do
             end if
         if(PR2.gt.1.0e-10) then
            if(abs(PR1/PR2).gt.(pi/2))then
             pr(j)=cos(pi-abs(PR1/PR2))
             else
             pr(j)=cos(PR1/PR2)
             endif
         else
           pr(j)=0.0
         end if
         end do
         mincc=maxval(pr)
         maxN=maxloc(pr)
         print*,maxN
        particle(i)%theta=ort60( 1,maxN(1))
        particle(i)%phi  =ort60( 2,maxN(1))
        particle(i)%omega=ort60( 3,maxN(1))
        particle(i)%PR   =mincc
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,   particle(i)%astigang,maxN

    end do
    deallocate(img2d,core_FTi,core_FTr,shift2d,a,b,a1,b1,a3)
    close(11)
    close(12)
    close(21)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,i5)
    return
end







