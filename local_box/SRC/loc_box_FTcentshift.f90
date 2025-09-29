


subroutine loc_box_FTcentshift
use common_icos_loc_box
use FFTW3SPACE
    implicit none
    integer i,ii,iii,i2,j,hx,ky,lz,imgID0,no0,no1,n1,n2,imgparticle0,ix,iy,ix0,iy0,ip1,ip2,Nparticle
    real theta,phi,omega,dr,rot(3,3),x0(3),x(3),x00(3)
    real matsym(3,3,60)
    real loc_ort(3,60),loc_cent(3,60),loc_cent0(3,60),rot0(3,3)
    character(len=54)::stckfile

    real,allocatable::loc_a(:),loc_b(:),loc_img(:,:)
    real,allocatable::a0(:),b0(:),a(:),b(:),shift2d(:),a00(:)
    real,allocatable::img2dstck(:)

    if(proc2d_icos_imgstck.eq.'y') then
        open(11,file=trim(icos_imgstck),form='unformatted',access='stream',status='old')
        read(11) mrc
    end if
    open(21,file=trim(local_imgstck), form='unformatted',access='stream')
    open(22,file=trim(local_ort))
!    open(23,file='test.mrc',form='unformatted',access='stream')
    allocate(a00(FFTsize**2),a0(FFTsize**2),b0(FFTsize**2),a(FFTsize**2),b(FFTsize**2),loc_img(sub_FFTsize**2,60),shift2d(FFTsize**2))

    call loc_maphead
    mrc%nz=60*totalparticle
    write(21) mrc

    imgID0=-100
    imgparticle0=-100
    matsym=.0
    call icos2f_sym_matrix(matsym,60)

    no0=-1
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))

    shift2d=-1
    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        if(mod((hx+ky),2).eq.0) shift2d(i2)=1
!        shift2d(i2)=(-1)**(hx+ky)
    end do

    x0(1)=loc_x;x0(2)=loc_y;x0(3)=loc_z
    Nparticle=0
    do i=1,totalparticle
        print*,i
        ! read raw particle image
        if(proc2d_icos_imgstck.eq.'y') then
            read(11) a0
        else

        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(11) mrc
            allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
            read(11) img2dstck
            close(11)
        end if
        a0(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        end if
        call img2d_norm(a0,FFTsize)
        a00=a0

        !calculate 60 equal orientations
        loc_ort=.0  
        call icosequalort(particle(i)%theta,particle(i)%phi,particle(i)%omega,loc_ort,60)

        !calculate 60 sub_images
        loc_cent=.0
        b0=0
        loc_img=.0
        do ii=1,60
            call Eularmatrix(loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii),rot)
            loc_cent(:,ii)=matmul(transpose(rot),x0)  !sub_img center in particle
            Nparticle=Nparticle+1
            write(22,900) particle(i)%stckfile,Nparticle,loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii), &
                   .0,.0,particle(i)%PR, particle(i)%imgID,ii,particle(i)%df1,particle(i)%df2,particle(i)%astigang,loc_cent(3,ii)*apix
            a=a0
            call img2d_centshift(a,FFTsize,-loc_cent(1,ii)-particle(i)%x,-loc_cent(2,ii)-particle(i)%y) !shift the sub_img cent to img center
            do iii=1,sub_FFTsize**2
                ky=(iii-1)/sub_FFTsize
                hx=mod(iii-1,sub_FFTsize)
                hx=hx-(sub_FFTsize-FFTsize)/2
                ky=ky-(sub_FFTsize-FFTsize)/2
                i2=ky*FFTsize+hx+1
                if((i2.ge.1).and.(i2.le.FFTsize**2)) loc_img(iii,ii)=a(i2)
            end do

        end do
        write(21) loc_img
        
    end do
    return
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f12.4)

end

 