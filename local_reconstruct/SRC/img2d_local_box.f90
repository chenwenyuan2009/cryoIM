subroutine img2d_local_box(img2d,FFTsize0,ix,iy,local_img2d,FFTsize)
implicit none
    integer FFTsize0,FFTsize,ix,iy,hx,ky,i,i2
    real img2d(FFTsize0**2),local_img2d(FFTsize**2)

    do i=1,FFTsize**2
        ky=(i-1)/FFTsize
        hx=mod(i-1,FFTsize)
        hx=hx+ix-FFTsize/2
        ky=ky+iy-FFTsize/2
        if(hx.lt.0) hx=hx+FFTsize0
        if(hx.ge.FFTsize0) hx=hx-FFTsize0
        if(ky.lt.0) ky=ky+FFTsize0
        if(ky.ge.FFTsize0) ky=ky-FFTsize0
        i2=ky*FFTsize0+hx+1
        if((i2.ge.1).and.(i2.le.FFTsize0**2)) local_img2d(i)=img2d(i2)
    end do
    return
end
 