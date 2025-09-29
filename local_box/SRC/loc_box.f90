! local boxer
subroutine loc_box
use common_icos_loc_box
use FFTW3SPACE
    implicit none
    integer i,ii,iii,i2,hx,ky,hx0,ky0,Nparticle
    real rot(3,3),x0(3)
    character(len=64)::stckfile
    real,allocatable::loc_ort(:,:),loc_cent(:,:),a0(:),a1(:),loc_img(:,:),img2dstck(:)

    allocate(loc_ort(3,Nsym),loc_cent(3,Nsym))
    allocate(a0(FFTsize**2),a1(sub_FFTsize**2),loc_img(sub_FFTsize**2,Nsym))

    if(proc2d_icos_imgstck.eq.'y') then
        open(11,file=trim(icos_imgstck),form='unformatted',access='stream',status='old')
        read(11) mrc
        do i=1,first-1
            call fseek(11,4*FFTsize**2,1)
        end do
    end if
    open(21,file=trim(local_imgstck), form='unformatted',access='stream')
    open(22,file=trim(local_ort))
    call loc_maphead
    mrc%nz=Nsym*totalparticle
    write(21) mrc

    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))

    x0(1)=loc_x;x0(2)=loc_y;x0(3)=loc_z
    Nparticle=0
    do i=first,last !1,totalparticle
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

        if(particle(i)%status.eq.'g') then
        call img2d_centshift(a0,FFTsize,-particle(i)%x,-particle(i)%y)
        call img2d_norm(a0,FFTsize)
        !calculate equal orientations
        loc_ort=.0  
        call equalort(particle(i)%theta,particle(i)%phi,particle(i)%omega,sym_rot,loc_ort,Nsym)

        !calculate sub_images
        loc_cent=.0
        loc_img=.0
        do ii=1,Nsym
            call Eularmatrix(loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii),rot)
            loc_cent(:,ii)=matmul(transpose(rot),x0)  !sub_img center in particle
            Nparticle=Nparticle+1
            hx0=int(loc_cent(1,ii))
            loc_cent(1,ii)=loc_cent(1,ii)-hx0
            ky0=int(loc_cent(2,ii))
            loc_cent(2,ii)=loc_cent(2,ii)-ky0
            hx0=hx0+FFTsize/2-sub_FFTsize/2
            ky0=ky0+FFTsize/2-sub_FFTsize/2
            do iii=1,sub_FFTsize**2
                ky=(iii-1)/sub_FFTsize+ky0
                hx=mod(iii-1,sub_FFTsize)+hx0
!                if((hx.gt.0).and.(hx.lt.FFTsize).and.(ky.gt.0).and.(ky.lt.FFTsize)) then
!                    i2=ky*FFTsize+hx+1
!                    loc_img(iii,ii)=a0(i2)
!                end if
                if(hx.lt.0) hx=hx+FFTsize
                if(hx.ge.FFTsize) hx=hx-FFTsize
                if(ky.lt.0) ky=ky+FFTsize
                if(ky.ge.FFTsize) ky=ky-FFTsize
                i2=ky*FFTsize+hx+1
                loc_img(iii,ii)=a0(i2)
            end do
            a1=loc_img(:,ii)
            call img2d_norm(a1,sub_FFTsize)
            loc_img(:,ii)=a1
            write(22,900) particle(i)%stckfile,Nparticle,loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii), &
                   loc_cent(1,ii),loc_cent(2,ii),particle(i)%PR, particle(i)%imgID,ii,  &
                   particle(i)%df1,particle(i)%df2,particle(i)%astigang,loc_cent(3,ii)*apix
        end do
        write(21) loc_img
        end if
    end do
    close(11)
    close(21)
    close(22)
    return
900 format(a64,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f12.4)

end

 