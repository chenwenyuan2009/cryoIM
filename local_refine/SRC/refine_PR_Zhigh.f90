

subroutine refine_PR_Zhigh
use common_loc_refine
    implicit none
!    integer pID,Nmax(1),item,Icycle,hx,ky,lz,imgID0,maxN(1),i1,i3,i2,nr,imgparticleID0,no1,no0
!    real cc,theta0,phi0,omega0,x0,y0,z0,maxCC,dangle,dtran,dr,PR,r2,df1,df2
!    real ort13(5,0:10),cc13(0:10)
!    character(len=64)::stckfile
!    real,allocatable:: a1(:),a2(:),a20(:),model2d(:,:,:),ctf2d(:),ctf2d_13(:,:),img3dr(:),img3di(:),shift3d(:),shift2d(:)
!    real,allocatable:: img2d(:),img2dstck(:)
!    real,allocatable:: b1(:),b2(:),b20(:),w_lp(:)
!
!    allocate(img2d(FFTsize0**2))
!    allocate(a1(FFTsize**2),a2(FFTsize**2),a20(FFTsize**2),w_lp(FFTsize**2))
!    allocate(b1(FFTsize**2),b2(FFTsize**2),b20(FFTsize**2),shift2d(FFTsize**2))
!    allocate(model2d(FFTsize**2,2,0:10))
!    allocate(img3dr(FFTsize**3),img3di(FFTsize**3),shift3d(FFTsize**3))
!    allocate(ctf2d(FFTsize**2),ctf2d_13(FFTsize**2,0:10))
    integer pID,Nmax(1),item,Icycle,hx,ky,lz,imgID0,maxN(1),i3,i2,nr,imgparticleID0,no0,no1,i
    real cc,theta0,phi0,omega0,x0,y0,maxCC,dangle,dtran,dr,PR,r2,df1,df2,theta1,phi1,omega1
    real ort11(5,0:10),cc11(0:10),z0
    character(len=64)::stckfile
    real,allocatable:: a1(:),a2(:),a20(:),model2d(:,:,:),ctf2d(:),img3dr(:),img3di(:),shift3d(:),shift2d(:)
    real,allocatable:: img2d(:),img2dstck(:),ortN(:,:),a12(:)
    real,allocatable:: b1(:),b2(:),b20(:),w_lp(:),b12(:),pr11(:)

    allocate(img2d(FFTsize0**2))
    allocate(a1(FFTsize**2),a2(FFTsize**2),a20(FFTsize**2),w_lp(FFTsize**2))
    allocate(b1(FFTsize**2),b2(FFTsize**2),b20(FFTsize**2),shift2d(FFTsize**2))
    allocate(model2d(FFTsize**2,2,0:10),a12(FFTsize**2),b12(FFTsize**2))
    allocate(img3dr(FFTsize**3),img3di(FFTsize**3),shift3d(FFTsize**3))
    allocate(ctf2d(FFTsize**2),ortN(3,Nsym),pr11(Nsym))

    open(32,file=trim(model3d),form='unformatted',access='stream',status='old')
    open(41,file=trim(local_newort))

    read(32) mrc
    read(32) img3dr  !read 3D model
    call img3d_softmask(img3dr,FFTsize,imgmask)
    close(32)

    shift3d=1.
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) shift3d(i3)=-1.
    end do
    shift2d=1.
    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        if(mod(hx+ky,2).ne.0) shift2d(i2)=-1.
    end do

    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        r2=float((hx-FFTsize/2)**2+(ky-FFTsize/2)**2)/2.0/maxRes**2
        w_lp(i2)=exp(-r2)
    end do

    img3dr=img3dr*shift3d
    img3di=.0
    call img3d_FFT(img3dr,img3di,FFTsize,1)
    img3dr=img3dr*shift3d
    img3di=img3di*shift3d
    no0=-1
    ctf2d=1.0
    stckfile='abcdef'
    imgID0=-100
    imgparticleID0=-100
    allocate(img2dstck(FFTsize0**2))
    if(proc2d%imgstck.eq.'y') then
        open(31,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(31) mrc
        do pID=1,first-1
            call fseek(11,4*FFTsize0**2,1)
        end do
    end if

    do pID=first,last
            no1=int(float((pID-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,5).eq.0)) then
            no0=no1
            write(*,*) no0,'%'
        end if
        if(proc2d%imgstck.eq.'y') then
            read(31) img2d
        else
            if(trim(stckfile).ne.trim(particle(pID)%stckfile)) then
                deallocate(img2dstck)
                stckfile=trim(particle(pID)%stckfile)
                open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
                read(11) mrc
                allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
                read(11) img2dstck
                close(11)
            end if
            if((imgparticleID0.ne.particle(pID)%imgparticleID).or.(particle(pID)%imgID.ne.imgID0)) then
                img2d(:)=img2dstck((particle(pID)%imgparticleID-1)*FFTsize0**2+1:particle(pID)%imgparticleID*FFTsize0**2)
                call img2d_norm(img2d,FFTsize0)
                imgparticleID0=particle(pID)%imgparticleID
            end if
        end if

        call img2d_local_box(img2d,FFTsize0,particle(pID)%local_hx0,particle(pID)%local_ky0,a20,FFTsize)
        call img2d_mask(a20,FFTsize,imgmask,particle(pID)%local_x0,particle(pID)%local_y0)
        !    call img2d_norm(a,FFTsize)

        a20=a20*shift2d
        b20=.0
        call img2d_FFT(a20,b20,FFTsize,FFTsize,1)
        a20=a20*w_lp
        b20=b20*w_lp

        theta0=particle(pID)%local_theta0
        phi0  =particle(pID)%local_phi0
        omega0=particle(pID)%local_omega0
        x0    =particle(pID)%local_x0
        y0    =particle(pID)%local_y0
        z0    =particle(pID)%local_z0
         if(trim(sign).eq.'+')then    
              df1=particle(pID)%df1+z0
              df2=particle(pID)%df2+z0
        elseif(trim(sign).eq.'-')then
              df1=particle(pID)%df1-z0
              df2=particle(pID)%df2-z0
        endif
        call getCTF2d(FFTsize,df1,df2,particle(pID)%astigang,apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)

!        do Icycle=1,Ncycle
!            dangle =refine_angle*(convF)**(Icycle-1)
!            dtran  =refine_shift*(convF)**(Icycle-1)
!            call ort_cent_11(theta0,phi0,omega0,x0,y0,dangle,dtran,ort13)
!            call cal_template_FT(img3dr,img3di,ort13,FFTsize,model2d)
!            do item=0,10
!                a1=model2d(:,1,item)
!                b1=model2d(:,2,item)
!                a2=a20
!                b2=b20
!                a1=a1*shift2d*CTF2d*w_lp
!                b1=b1*shift2d*CTF2d*w_lp
!                minRes=-1.0
!                call cal_PR_FT(a1,b1,a2,b2,FFTsize,minRes,maxRes,PR)
!                cc13(item)=PR
!            end do
!            maxCC=minval(cc13)
!            maxN =minloc(cc13)
!            theta0=ort13(1,maxN(1)-1)
!            phi0  =ort13(2,maxN(1)-1)
!            omega0=ort13(3,maxN(1)-1)
!            x0    =ort13(4,maxN(1)-1)
!            y0    =ort13(5,maxN(1)-1)
!            
!        end do
     if(sym.ne.'c1')then
        theta1=particle(pID)%local_theta0
        phi1  =particle(pID)%local_phi0
        omega1=particle(pID)%local_omega0
        x0    =particle(pID)%local_x0
        y0    =particle(pID)%local_y0
        if(sym.eq.'asym'.or.sym.eq.'ASYM')then
        phi1  =0.0
        call icos_equalortasym(theta1,phi1,omega1,ortN)
        else
        call icos_equalort(theta1,phi1,omega1,ortN)
        endif
        do i=1,Nsym
           call FFT3d_slice(img3dr,img3di,FFTsize,ortN(1,i),ortN(2,i),ortN(3,i),a12,b12)
           call img2d_FTphase(a12,b12,FFTsize,x0,y0)
                a12=a12*shift2d*CTF2d*w_lp
                b12=b12*shift2d*CTF2d*w_lp
                a2=a20
                b2=b20
                minRes=-1.0
                call cal_PR_FT(a12,b12,a2,b2,FFTsize,minRes,maxRes,PR)
                pr11(i)=PR
         end do
            maxCC=minval(pr11)
            maxN =minloc(pr11)
            theta0=ortN(1,maxN(1))
            phi0  =ortN(2,maxN(1))
            omega0=ortN(3,maxN(1))
        if(Ncycle.ne. 0)then
        do Icycle=1,Ncycle
            dangle =refine_angle*(convF)**(Icycle-1)
            dtran  =refine_shift*(convF)**(Icycle-1)
            call ort_cent_11(theta0,phi0,omega0,x0,y0,dangle,dtran,ort11)
            call cal_template_FT(img3dr,img3di,ort11,FFTsize,model2d)
            do item=0,10
                a1=model2d(:,1,item)
                b1=model2d(:,2,item)
                a2=a20
                b2=b20
                a1=a1*shift2d*CTF2d*w_lp
                b1=b1*shift2d*CTF2d*w_lp
                minRes=-1.0
                call cal_PR_FT(a1,b1,a2,b2,FFTsize,minRes,maxRes,PR)
                cc11(item)=PR
            end do
            maxCC=minval(cc11)
            maxN =minloc(cc11)
            theta0=ort11(1,maxN(1)-1)
            phi0  =ort11(2,maxN(1)-1)
            omega0=ort11(3,maxN(1)-1)
            x0    =ort11(4,maxN(1)-1)
            y0    =ort11(5,maxN(1)-1)
        end do
        endif 
        else
        theta0=particle(pID)%local_theta0
        phi0  =particle(pID)%local_phi0
        omega0=particle(pID)%local_omega0
        x0    =particle(pID)%local_x0
        y0    =particle(pID)%local_y0
        do Icycle=1,Ncycle
            dangle =refine_angle*(convF)**(Icycle-1)
            dtran  =refine_shift*(convF)**(Icycle-1)
            call ort_cent_11(theta0,phi0,omega0,x0,y0,dangle,dtran,ort11)
            call cal_template_FT(img3dr,img3di,ort11,FFTsize,model2d)
            do item=0,10
                a1=model2d(:,1,item)
                b1=model2d(:,2,item)
                a2=a20
                b2=b20
                a1=a1*shift2d*CTF2d*w_lp
                b1=b1*shift2d*CTF2d*w_lp
                minRes=-1.0
                call cal_PR_FT(a1,b1,a2,b2,FFTsize,minRes,maxRes,PR)
                cc11(item)=PR
            end do
            maxCC=minval(cc11)
            maxN =minloc(cc11)
            theta0=ort11(1,maxN(1)-1)
            phi0  =ort11(2,maxN(1)-1)
            omega0=ort11(3,maxN(1)-1)
            x0    =ort11(4,maxN(1)-1)
            y0    =ort11(5,maxN(1)-1)
        end do
        endif 

        particle(pID)%local_dtheta=theta0-particle(pID)%local_theta0
        particle(pID)%local_dphi  =phi0  -particle(pID)%local_phi0
        particle(pID)%local_domega=omega0-particle(pID)%local_omega0
        particle(pID)%local_dx    =x0-particle(pID)%local_x0
        particle(pID)%local_dy    =y0-particle(pID)%local_y0
        particle(pID)%local_dz    =z0-particle(pID)%local_z0
        particle(pID)%local_theta =theta0
        particle(pID)%local_phi   =phi0
        particle(pID)%local_omega =omega0
        particle(pID)%local_x     =x0
        particle(pID)%local_y     =y0
        particle(pID)%local_z     =z0
        particle(pID)%local_PR    =maxCC
        particle(pID)%local_dPR   =particle(pID)%local_PR-particle(pID)%local_PR0

        write(41,901) particle(pID)%stckfile,particle(pID)%particleID,particle(pID)%theta,particle(pID)%phi,particle(pID)%omega, &
                particle(pID)%x,particle(pID)%y,particle(pID)%local_PR,particle(pID)%imgID,particle(pID)%imgparticleID,  &
                particle(pID)%df1,particle(pID)%df2,particle(pID)%astigang,  &
                particle(pID)%local_imgID,particle(pID)%local_hx0,particle(pID)%local_ky0, &
                particle(pID)%local_theta,particle(pID)%local_phi,particle(pID)%local_omega,  &
                particle(pID)%local_x,particle(pID)%local_y,particle(pID)%local_z,particle(pID)%local_PR,  &
                particle(pID)%local_dtheta,particle(pID)%local_dphi,particle(pID)%local_domega,  &
                particle(pID)%local_dx,particle(pID)%local_dy,particle(pID)%local_dz,particle(pID)%local_dPR


    end do
    close(31)
    close(41)

901 format(a64,I9,F10.4,F10.4,F10.4,  &
           F10.4,F10.4,F10.4,I11,I7,  &
           F15.4,F15.4,F10.4, &
           i9,i6,i6,  &
           f10.4,f10.4,f10.4,  &
           f10.4,f10.4,f11.4,f10.4,  &
           f10.4,f10.4,f10.4,  &
           f10.4,f10.4,f10.4,f10.4)
!    deallocate(a1,a2,a20,model2d,ctf2d,img3dr,img3di,shift3d,shift2d)
!    deallocate(b1,b2,b20)
    deallocate(a1,a2,a20,model2d,ctf2d,img3dr,img3di,shift3d,shift2d)
    deallocate(b1,b2,b20)
    return
 end subroutine
