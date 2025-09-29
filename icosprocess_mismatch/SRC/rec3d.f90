subroutine rec3d
   use common_icosproc
   USE IFPORT
   use FFTW3SPACE
   implicit none
    integer pid,i2,hx,ky,lz,n1,n2,n3,i3,i4,imgID0,no0,no1
    integer nx0,ny0,nz0,nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz
    integer ix,iy,iz,usedparticle
    integer arrayID2,arrayID3,isym
    real theta,phi,omega,dr,di,rot(3,3),x0(3),x(3),x00(3),boundX
    real matsym(3,3,60),defocus,RI,xbit,ybit,zbit,addweight,weight
    real symoperate(3,3),calculatesinc,getXYCTF,ctf0,Rseed

    real(kind=4),allocatable::wg(:),wg1(:)
    real,allocatable::a(:),ctf2d(:)
    complex,allocatable::FFT2d(:),FFT3d(:),FFT3d1(:),A3d(:)

    allocate(FFT3d(FFTsize**3))
    allocate(wg(FFTsize**3))
    allocate(a(FFTsize**2))
    allocate(ctf2d(FFTsize**2))


    open(11,file=trim(imgstck),form='binary',status='old')
    n1=FFTsize
    n2=FFTsize
    n3=FFTsize

    boundX=0.3*FFTsize/2

    FFT3d=cmplx(0.0,0.0)
    Rseed=100.0
    wg=1.0/SSNR
    imgID0=-100
    ctf2d=1.0
    usedparticle=0
    matsym=.0
    defocus=0.0
    matsym(1,1,1)=1;matsym(2,2,1)=1;matsym(3,3,1)=1

!    if(sym.eq.'T') call sym_T_222(matsym)
!    if(sym.eq.'D') call sym_D2(matsym)
!    if(sym.eq.'O') call sym_octa_422(matsym)
!    if(sym.eq.'I') call sym_icos_222(matsym)


    read(11) mrc

    offset=4*FFTsize**2
    ipos=1

    do pid=1,first-1
        istat=fseek(11,offset,ipos)
    end do

    allocate(FFT2d(FFTsize**2))

    call intialplanCmplx2d(FFT2d,FFT2d,n1,n2)
    call intialplanCmplx3d(FFT3d,FFT3d,n1,n2,n3)

    mrc%nz=wholegroup

    do pid=first,last

      if(particle(pid)%good) then
        read(11) a
        usedparticle=usedparticle+1

         write(*,FMT="(a,1x,T10,F6.2,'%',$)") char(13),real(usedparticle)/real(totalparticle)*100.0

        if(needimgmask.eq.'y') call img2d_mask(a,particle(pid)%centx,particle(pid)%centy)
        if(neednorm.eq.'y') call img2d_norm(a)

        if((NeedCtf.eq.'y').and.(particle(pid)%group.ne.imgID0)) then
            if(particle(pid)%df2.gt. 0.001) then
                defocus=0.5*particle(pid)%df1+0.5*particle(pid)%df2
            else
                defocus=particle(pid)%df1
            end if
            call getCTF2d(ctf2d,defocus)
            imgID0=particle(pid)%group
        end if

        FFT2d=cmplx(a,0.0)

        call shift2dComplex(FFT2d,n1,n2)
        call FFTcmplx2d_f(FFT2d,FFT2d,n1,n2)

         FFT2d=conjg(FFT2d)

        call FTcentphaseshift(FFT2d,particle(pid)%centx,particle(pid)%centy)



        theta=particle(pid)%theta
        phi  =-particle(pid)%phi
        omega=particle(pid)%omega+180.0
        call oular_mat(rot,theta,phi,omega)
        ctf0=-1.0


        do ky=-maxFR,maxFR-1
            do hx=-maxFR,maxFR-1
                dr=sqrt(float(hx**2+ky**2))
                if(dr.le.maxFR) then
                i2=arrayID2(hx+FFTsize/2,ky+FFTsize/2)
                if(NeedCtf.eq.'y') ctf0=ctf2d(i2)
                    x(1)=hx
                    x(2)=ky
                    x(3)=.0
                    x0=matmul(rot,x)
!************************************************************************************
                    do isym=1,nsym

                       symoperate=matsym(:,:,isym)
                       x=matmul(symoperate,x0)
!************************************************************************************
                    if(sqrt(X(1)**2.0+X(2)**2.0+X(3)**2.0) .ge. maxFR) cycle

                    nx0=floor(x(1));ny0=floor(x(2));nz0=floor(x(3))
                    nx1=nx0;ny1=ny0;nz1=nz0
                    nx2=nx0+1;ny2=ny0+1;nz2=nz0+1

                  do iz=nz1,nz2         !
                        do iy=ny1,ny2
                            do ix=nx1,nx2
                                 xbit=abs(ix-X(1))
                                 ybit=abs(iy-X(2))
                                 zbit=abs(iz-X(3))
                                 di=addweight(xbit,ybit,zbit)
                                 i3=arrayID3(ix+FFTsize/2,iy+FFTsize/2,iz+FFTsize/2)
                                 weight=di*ctf0
                                 FFT3d(i3)=FFT3d(i3)+(FFT2d(i2))*weight
                                 wg(i3)=wg(i3)+weight*ctf0
                            end do
                        end do
                    end do
                    end do
                end if
              end do
           end do
     else
        istat=fseek(11,offset,ipos)
     end if
    end do


    FFT3d=FFT3d/wg


    deallocate(FFT2d)
    write(*,*)
    print*,'perform IFFT ...'

    FFT3d=conjg(FFT3d)
    call shift3dComplex(FFT3d,n1,n2,n3)
    call FFTcmplx3d_B(FFT3d,FFT3d,n1,n2,n3)
    call shift3dComplex(FFT3d,n1,n2,n3)
    call maphead
    wg=real(FFT3d)
    call img3d_mask(wg)
    call img3d_norm(wg)


    open(21,file=trim(result3d),form='binary')
    write(21) mrc
    write(21) wg
    close(11)
    close(21)

    deallocate(wg,FFT3d,ctf2d)

    return
end subroutine

