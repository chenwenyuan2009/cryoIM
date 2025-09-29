
!cosine file in real space
subroutine img3d_softmask(img3d,FFTsize,mask)
    implicit none
    integer FFTsize,i3,ix,iy,iz,halfw,num
    real pi,r,mask,innR,outR,summ,edge
    real img3d(FFTsize**3)

    halfw=6
    pi=4.0*atan2(1.0,1.0)
    innR=mask-halfw/2
    if(innR.le.0) innR=0
    outR=mask+halfw/2
    if(outR.gt.FFTsize/2) outR=FFTsize/2

    num=0
    summ=0
    do i3=1,FFTsize**3
        iz=(i3-1)/FFTsize**2-FFTsize/2
        if(iz.le.outR) then
            iy=mod(i3-1,FFTsize**2)/FFTsize-FFTsize/2
            ix=mod(mod(i3-1,FFTsize**2),FFTsize)-FFTsize/2
            r=sqrt(float(ix**2+iy**2+iz**2))
            if((r.gt.innR).and.(r.le.outR)) then
                summ=img3d(i3)+summ
                num=num+1
            end if
        end if
    end do
    if(num.ne.0) summ=summ/float(num)
    print*,'soft mask:',summ,num

    do i3=1,FFTsize**3
        iz=(i3-1)/FFTsize**2-FFTsize/2
        iy=mod(i3-1,FFTsize**2)/FFTsize-FFTsize/2
        ix=mod(mod(i3-1,FFTsize**2),FFTsize)-FFTsize/2
        r=sqrt(float(ix**2+iy**2+iz**2))
        if((r.gt.innR).and.(r.le.outR)) then
            edge=(1.0+cos(pi*(r-innR)/float(halfw)))/2.0
            img3d(i3)=img3d(i3)*edge+(1.0-edge)*summ
         end if
         if(r.gt.outR) img3d(i3)=summ
     end do

    return
end subroutine
