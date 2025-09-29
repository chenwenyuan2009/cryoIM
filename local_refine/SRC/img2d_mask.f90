!cosine file in real space
subroutine img2d_mask(Density,nx,mask,OrigX,OrigY)
    implicit none
    integer i,iy,ix,num,nx,halfw
    real mask,pi
    real R1,R2,Ra,xx,yy,edge,OrigX,OrigY,x0,y0,MaskR0
    double precision summ
    real Density(nx*nx)

    MaskR0=mask
    halfw=6
    pi=4.0*atan2(1.0,1.0)
    R1=MaskR0-halfw/2
    if(R1.le.0) R1=0
    R2=MaskR0+HalfW/2
    if(R2.gt.nx/2) R2=nx/2

    summ=.0
    num=0
    x0=OrigX+nx/2
    y0=OrigY+nx/2
    do i=1,nx*nx
        iy=(i-1)/nx
        ix=mod(i-1,nx)
        yy=(-y0+iy)**2
        xx=(-x0+ix)**2
        ra=sqrt(yy+xx)
        if((Ra.gt.R1).and.(Ra.le.R2)) then
            summ=Density(i)+summ
            num=num+1
        end if
    end do
    if(num.ne.0) summ=summ/float(num)

    !mask the 2D image according to cosine edge.
    do i=1,nx*nx
        iy=(i-1)/nx
        ix=mod(i-1,nx)
        yy=(-y0+iy)**2
        xx=(-x0+ix)**2
        ra=sqrt(yy+xx)
        if((Ra.gt.R1).and.(Ra.le.R2)) then
            edge=(1.0+cos(pi*(Ra-R1)/float(halfw)))/2.0
            Density(i)=Density(i)*edge+(1.0-edge)*summ
         end if
         if(Ra.gt.R2) Density(i)=summ
     end do
     return
end subroutine




