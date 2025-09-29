subroutine mapmask
    use common_icosproc
    implicit none
    integer i3,i,j,hx,ky,lz,z,i2,x11,y11,z11
    real r
    real,allocatable::a3(:)
    real,allocatable::core_FTr(:),core_FTi(:)
    allocate(core_FTr(FFTsize**3),a3(FFTsize**2),core_FTi((FFTsize/2)**3))

    print*,'reading 3D map...'
    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    open(12,file=trim(result3d),form='unformatted',access='stream')
    read(11) mrc
    print*,FFTsize
    print*,clipsize
    do i=1,FFTsize
      read(11) a3
      do j=1,FFTsize**2
        core_FTr((i-1)*FFTsize**2+j)=a3(j)
      end do
    end do
    close(11)
!    if(proc2d%clipsize.eq.'y')then
    mrc%nx=FFTsize/2
    mrc%ny=FFTsize/2
    mrc%nz=FFTsize/2
    write(12) mrc
!    else
!    write(12) mrc
!    endif
    !read(*,*) x11,y11,z11

!    endif
!    do lz=z11+FFTsize/2-clipsize/2,z11+FFTsize/2+clipsize/2-1
!        do ky=y11+FFTsize/2-clipsize/2,y11+FFTsize/2+clipsize/2-1
!            do hx=x11+FFTsize/2-clipsize/2,x11+FFTsize/2+clipsize/2-1
!               !r=sqrt(float(lz-FFTsize/2-1)**2+float(ky-FFTsize/2-1)**2+float(hx-FFTsize/2-1)**2)
!               i3=(lz-1)*FFTsize**2+(ky-1)*FFTsize+hx
!               i2=(lz+clipsize/2-FFTsize/2-z11)*(clipsize/2)**2+(ky-y11-FFTsize/2+clipsize/2-1)*(clipsize/2)&
!                 +hx-FFTsize/2+clipsize/2-x11+1
!!               if(r.lt.capsid_innerR.or.r.ge.capsid_outerR) then
!               core_FTi(i2)=core_FTr(i3)
!!               endif
!            end do
!        end do
!    end do
    do lz=1,FFTsize/2
        do ky=1,FFTsize/2
            do hx=1,FFTsize/2
               !r=sqrt(float(lz-FFTsize/2-1)**2+float(ky-FFTsize/2-1)**2+float(hx-FFTsize/2-1)**2)
               i3=(lz-1)*FFTsize**2+(ky-1)*FFTsize+hx
               i2=(lz-1)*(FFTsize/2)**2+(ky-1)*FFTsize/2+hx
!               if(r.lt.capsid_innerR.or.r.ge.capsid_outerR) then
               core_FTi(i2)=core_FTr(i3)
!               endif
            end do
        end do
    end do
!    if(proc2d%clipsize.eq.'y')then
!    do z=1,clipsize
!        do j=1,clipsize
!           do i=1,clipsize
!            i3=(z+(FFTsize-clipsize)/2)*FFTsize**2+(j+(FFTsize-clipsize)/2)*FFTsize+i+(FFTsize-clipsize)/2
!            i2=(z-1)*clipsize**2+(j-1)*clipsize+i
!            core_FTi(i2)=core_FTr(i3)
!            enddo
!        end do
!    end do
    write(12) core_FTi
!    else
!    write(12) core_FTr
!    endif
    close(12)
    return
end

