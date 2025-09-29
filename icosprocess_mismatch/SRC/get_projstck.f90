subroutine get_projstck
    use common_icosproc
    real,allocatable:: imgproj(:)

    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    open(21,file='proj.stck',form='unformatted',access='stream')

    read(11)mrc
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

