subroutine genomeimgstck
use common_icosproc
    implicit none
    integer i3,i,hx,ky,lz,imgID0,no0,no1,j
    real,allocatable::a0(:),a1(:),a2(:),img2dstck(:),a(:),a3(:)
    real,allocatable::capsid_FTr(:),capsid_FTi(:),gene2d(:),coremodel_FTr(:),coremodel_FTi(:)

    allocate(a0(FFTsize**2),a1(FFTsize**2),a2(FFTsize**2),gene2d(geneFFTsize**2))
    allocate(capsid_FTr(FFTsize**3),capsid_FTi(FFTsize**3),a(FFTsize**2))
    allocate(coremodel_FTr(FFTsize**3),coremodel_FTi(FFTsize**3),a3(FFTsize**2))

    print*,'reading 3D map...'
    open(11,file=trim(model3d), form='unformatted',access='stream',status='old')
    open(12,file=trim(result3d),form='unformatted',access='stream',status='replace')
    open(31,file=trim(imgstck), form='unformatted',access='stream',status='old')
    if((proc3d%outputprojstck.eq.'y').or.(proc3d%outputprojstck.eq.'Y')) &
        open(13,file='proj.stck',form='unformatted',access='stream',status='replace')
    if((proc3d%outputrawimgstck.eq.'y').or.(proc3d%outputrawimgstck.eq.'Y')) &
        open(14,file='rawimg.stck',form='unformatted',access='stream',status='replace')
    open(21,file=trim(newortfile))
    if((proc3d%core3d.eq.'y').or.(proc3d%core3d.eq.'Y')) then
        open(33,file=trim(core3d),form='unformatted',access='stream',status='old')
        read(33) mrc
            do i=1,FFTsize
               read(33) a
                  do j=1,FFTsize**2
                     coremodel_FTr((i-1)*FFTsize**2+j)=a(j)
                  end do
            end do
        close(33)
        do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
          if(mod(hx+ky+lz,2).ne.0) then
            coremodel_FTr(i3)=-coremodel_FTr(i3)
          end if
        end do
        coremodel_FTi=.0

        print*,'preform 3D FFT...'
        call img3d_FFT(coremodel_FTr,coremodel_FTi,FFTsize,FFTsize,FFTsize,1)
        do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            coremodel_FTr(i3)=-coremodel_FTr(i3)
            coremodel_FTi(i3)=-coremodel_FTi(i3)
        end if
        end do
    endif

    read(11) mrc
    do i=1,FFTsize
      read(11) a
      do j=1,FFTsize**2
        capsid_FTr((i-1)*FFTsize**2+j)=a(j)
      end do
    end do
    close(11)
    mrc%nz=last-first+1
    if((proc3d%outputprojstck.eq.'y')  .or.(proc3d%outputprojstck.eq.'Y'))   write(13) mrc
    if((proc3d%outputrawimgstck.eq.'y').or.(proc3d%outputrawimgstck.eq.'Y')) write(14) mrc
    mrc%nx=geneFFTsize
    mrc%ny=geneFFTsize
    write(12) mrc

    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            capsid_FTr(i3)=-capsid_FTr(i3)
        end if
    end do
    capsid_FTi=.0

    print*,'preform 3D FFT...'
    call img3d_FFT(capsid_FTr,capsid_FTi,FFTsize,FFTsize,FFTsize,1)
    do i3=1,FFTsize**3
        lz=(i3-1)/FFTsize**2
        ky=mod(i3-1,FFTsize**2)/FFTsize
        hx=mod(mod(i3-1,FFTsize**2),FFTsize)
        if(mod(hx+ky+lz,2).ne.0) then
            capsid_FTr(i3)=-capsid_FTr(i3)
            capsid_FTi(i3)=-capsid_FTi(i3)
        end if
    end do

    print*,'Extracting gene images...'
    no0=-1
    read(31) mrc
    do i=1,first-1
        call fseek(31,4*FFTsize**2,1)
    end do
    do i=first,last
        print*,i
        read(31) a2

        call img2d_centshift_FT(a2,FFTsize,particle(i)%x,particle(i)%y)
        particle(i)%x=.0
        particle(i)%y=.0

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            if(particle(i)%imgID.ne.imgID0) then
                call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if
        if((proc3d%core3d.eq.'y').or.(proc3d%core3d.eq.'Y')) then
            call project2d_FT3dTo2d(coremodel_FTr,coremodel_FTi,FFTsize,particle(i)%theta,&
                 particle(i)%phi,particle(i)%omega,a3)
            call img2d_norm(a3,FFTsize)
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y'))call img2d_applyCTF(a3,FFTsize,ctf2d)
            call img2d_lowpass(a3,FFTsize,maxR)
            call img2d_norm(a3,FFTsize)

            call project2d_FT3dTo2d(capsid_FTr,capsid_FTi,FFTsize,particle(i)%theta,particle(i)%phi,particle(i)%omega,a1)
            call img2d_norm(a1,FFTsize)
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y'))call img2d_applyCTF(a1,FFTsize,ctf2d)

            call img2d_lowpass(a1,FFTsize,maxR)
            call img2d_lowpass(a2,FFTsize,maxR)
            call img2d_norm(a1,FFTsize)
            call img2d_norm(a2,FFTsize)
        if((proc3d%outputprojstck.eq.'y')  .or.(proc3d%outputprojstck.eq.'Y'))   write(13) a1
        if((proc3d%outputrawimgstck.eq.'y').or.(proc3d%outputrawimgstck.eq.'Y')) write(14) a2
        call img2d_subproj_new(a1,a2,a3,FFTsize,capsid_innerR,capsid_outerR)
        if(proc2d%geneFFTsize.eq.'y') then
            call img2d_clip(a3,FFTsize,geneFFTsize,gene2d)
        else
            gene2d=a3
        end if
        call img2d_norm(gene2d,geneFFTsize)
        write(12) gene2d
        goto 10000
        endif

        call project2d_FT3dTo2d(capsid_FTr,capsid_FTi,FFTsize,particle(i)%theta,particle(i)%phi,particle(i)%omega,a1)
        call img2d_norm(a1,FFTsize)
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) call img2d_applyCTF(a1,FFTsize,ctf2d)

        call img2d_lowpass(a1,FFTsize,maxR)
        call img2d_lowpass(a2,FFTsize,maxR)
        call img2d_norm(a1,FFTsize)
        call img2d_norm(a2,FFTsize)
        if((proc3d%outputprojstck.eq.'y')  .or.(proc3d%outputprojstck.eq.'Y'))   write(13) a1
        if((proc3d%outputrawimgstck.eq.'y').or.(proc3d%outputrawimgstck.eq.'Y')) write(14) a2
        call img2d_subproj(a1,a2,FFTsize,capsid_innerR,capsid_outerR)
        if(proc2d%geneFFTsize.eq.'y') then
            call img2d_clip(a1,FFTsize,geneFFTsize,gene2d)
        else
            gene2d=a1
        end if
        call img2d_norm(gene2d,geneFFTsize)
        write(12) gene2d
10000   particle(i)%PR=cos(particle(i)%PR*3.1415926/180.0)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
        particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1,       particle(i)%df2,   particle(i)%astigang

    end do
    deallocate(a,a0,a1,a2,a3)
    deallocate(gene2d,capsid_FTr,capsid_FTi,coremodel_FTi,coremodel_FTr)
    close(12)
    close(13)
    close(14)
    close(21)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4)
    return
end



