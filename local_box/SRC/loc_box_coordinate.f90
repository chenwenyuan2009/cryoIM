! local boxer
subroutine loc_box_coordinate
use common_icos_loc_box
use FFTW3SPACE
    implicit none
    integer i,ii,iii,i2,hx,ky,hx0,ky0,Nparticle
    real rot(3,3),x0(3)
    character(len=64)::stckfile
    real,allocatable::loc_ort(:,:),loc_cent(:,:)

    allocate(loc_ort(3,Nsym),loc_cent(3,Nsym))

    open(22,file=trim(local_ort))
    x0(1)=loc_x;x0(2)=loc_y;x0(3)=loc_z
    Nparticle=0
    do i=first,last !1,totalparticle
        print*,i
        if(particle(i)%status.eq.'g') then
            loc_ort=.0  
            call equalort(particle(i)%theta,particle(i)%phi,particle(i)%omega,sym_rot,loc_ort,Nsym)
            loc_cent=.0
            do ii=1,Nsym
                call Eularmatrix(loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii),rot)
                loc_cent(:,ii)=matmul(transpose(rot),x0)+FFTsize/2  !sub_img center in particle
                loc_cent(1,ii)=loc_cent(1,ii)+particle(i)%x
                loc_cent(2,ii)=loc_cent(2,ii)+particle(i)%y
                Nparticle=Nparticle+1
                hx0=int(loc_cent(1,ii))
                loc_cent(1,ii)=loc_cent(1,ii)-hx0
                ky0=int(loc_cent(2,ii))
                loc_cent(2,ii)=loc_cent(2,ii)-ky0
                loc_cent(3,ii)=(loc_cent(3,ii)-FFTsize/2)*apix
!                hx0=hx0+FFTsize/2-sub_FFTsize/2
!                ky0=ky0+FFTsize/2-sub_FFTsize/2
                write(22,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                   particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                   particle(i)%df1*10000.0,particle(i)%df2*10000.0,particle(i)%astigang,  &
                   ii,hx0,ky0,loc_ort(1,ii),loc_ort(2,ii),loc_ort(3,ii),loc_cent(1,ii),loc_cent(2,ii),loc_cent(3,ii),particle(i)%PR
            end do
        end if
    end do
    close(11)
    close(21)
    close(22)
    return
900 format(a64,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f15.4,f15.4,f10.4,i9,i6,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f11.4,f11.4)

end

 