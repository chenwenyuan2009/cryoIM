subroutine test_applyctf
    use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,imgID0,Nmax(1),no0,no1,NNmax(1),nproc
    integer ix0,iy0,ix,iy,n1,ip1,ip2,x0,y0,ix00,iy00,jj,ii
    real dx,dy,cc,dn,dm,RotMat(3,3),delts,dang,dtran,ortcent(5,11),dfcent(3,7),dfpr(7)
    real ort60(3,5),cent60(3,5),Cmax,z(-1:1,-1:1),initiald(3,1),newd(3,1),cent11(3,11)
    real dforigin(3)
    real,allocatable::a0(:),a1(:),a2(:),a(:),a3(:),a4(:)
    real,allocatable::core_FTr(:),core_FTi(:),img2dstck(:),loc_img2d2(:),loc_proj2d2(:)
    real,allocatable::loc_img2d(:,:),loc_proj2d(:,:),loc_img2d1(:,:),loc_proj2d1(:,:)
    real(kind=8) Pr,phaseRedisual

    allocate(a1(FFTsize**2))
    allocate(core_FTr(FFTsize**3),core_FTi(FFTsize**3),a(FFTsize**2))

     open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
     open(20,file=trim(result3d),form='unformatted',access='stream')
    read(11) mrc
    mrc%nz=last*2
    write(20)mrc
    do i=1,FFTsize
      read(11) a
      do j=1,FFTsize**2
        core_FTr((i-1)*FFTsize**2+j)=a(j)
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
    print*,'FFTsize',FFTsize
        imgID0=-1

        do i=first,last
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
        if(particle(i)%imgID.ne.imgID0) then
        call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,apix, &
                              CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                imgID0=particle(i)%imgID
            end if
        end if
        call project2d_FT3dTo2d(core_FTr,core_FTi,FFTsize,particle(i)%theta,particle(i)%phi&
            ,particle(i),a1)
        write(20) a1
        call img2d_applyCTF(a1,FFTsize,ctf2d)
        write(20) a1
        enddo

        deallocate(a1,core_FTi,core_FTr,a,ctf2d)
        end
