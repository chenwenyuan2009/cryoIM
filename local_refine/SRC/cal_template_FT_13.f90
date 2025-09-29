subroutine cal_template_FT_13(img3dr,img3di,ort11,nx,model2d)
implicit none
    integer i1,nx
    real theta,phi,omega,x,y
    real img3dr(nx**3),img3di(nx**3),model2d(nx**2,2,0:12),ort11(6,0:12)
    real,allocatable::FT2da(:),FT2db(:),a(:),b(:)
    allocate(FT2da(nx**2),FT2db(nx**2),a(nx**2),b(nx**2))

    theta=ort11(1,0)
    phi  =ort11(2,0)
    omega=ort11(3,0)
    call FFT3d_slice(img3dr,img3di,nx,theta,phi,omega,FT2da,FT2db)
!    model2d(:,1,0)=FT2da(:)
!    model2d(:,2,0)=FT2db(:)

    do i1=0,4
        a(:)=FT2da !model2d(:,1,0)
        b(:)=FT2db !model2d(:,2,0)
        call img2d_FTphase(a,b,nx,ort11(4,i1),ort11(5,i1))
        model2d(:,1,i1)=a(:)
        model2d(:,2,i1)=b(:)
    end do
    do i1=5,6
        model2d(:,1,i1)=FT2da
        model2d(:,2,i1)=FT2db
    end do

    do i1=7,12
        theta=ort11(1,i1)
        phi  =ort11(2,i1)
        omega=ort11(3,i1)    
        call FFT3d_slice(img3dr,img3di,nx,theta,phi,omega,FT2da,FT2db)
        call img2d_FTphase(FT2da,FT2db,nx,ort11(4,i1),ort11(5,i1))
        model2d(:,1,i1)=FT2da(:)
        model2d(:,2,i1)=FT2db(:)
    end do
    deallocate(FT2da,FT2db,a,b)
    return
   end

