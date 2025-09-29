
subroutine img3d_icosavg(img3d0,nx,r0)
    use common_icosproc
    implicit none
    integer nx,i3,hx,ky,lz,isym,ix,iy,iz,n1,n2,n3,n4,n5,n6,n7,n8,no0,no1
    integer i3b,i3c,i3d
    real img3d0(nx**3),symmat(3,3,60),matsymc12(3,3,12),matsymc5a(3,3,5),matsymc6(3,3,6),matsymc3(3,3,3)
    real r,r0,x0(3),x(3),ax,ay,az,matsymc15(3,3,15),matsymc4(3,3,4),matsymc8(3,3,8),matsymc2(3,3,2)
    real,allocatable::img3d1(:)
    allocate(img3d1(nx**3))

    print*,'performing real space icos average...'
!    symmat=.0
!    nsym  =1
    img3d1=.0
!    nsym  =60
!    call icos2f_sym_mat(symmat)

    matsymc12=.0
    matsymc12(1,1,1)=1;matsymc12(2,2,1)=1;matsymc12(3,3,1)=1
    matsymc15=.0
    matsymc15(1,1,1)=1;matsymc15(2,2,1)=1;matsymc15(3,3,1)=1
    matsymc6=.0
    matsymc6(1,1,1)=1;matsymc6(2,2,1)=1;matsymc6(3,3,1)=1
    matsymc3=.0
    matsymc3(1,1,1)=1;matsymc3(2,2,1)=1;matsymc3(3,3,1)=1
    matsymc4=.0
    matsymc4(1,1,1)=1;matsymc4(2,2,1)=1;matsymc4(3,3,1)=1
    matsymc5a=.0
    matsymc5a(1,1,1)=1;matsymc5a(2,2,1)=1;matsymc5a(3,3,1)=1
    matsymc2=.0
    matsymc2(1,1,1)=1;matsymc2(2,2,1)=1;matsymc2(3,3,1)=1
    matsymc8=.0
    matsymc8(1,1,1)=1;matsymc8(2,2,1)=1;matsymc8(3,3,1)=1
    if(sym.eq.'c12')  call icos2f_sym_matc12(matsymc12)
    if(sym.eq.'c15')   call icos2f_sym_matc15(matsymc15)
    if(sym.eq.'c5a')   call icos2f_sym_matc5a(matsymc5a)
    if(sym.eq.'c6')   call icos2f_sym_matc6(matsymc6)
    if(sym.eq.'c3') call icos2f_sym_matc3(matsymc3)
    if(sym.eq.'c4') call icos2f_sym_matc4(matsymc4)
    if(sym.eq.'c8') call icos2f_sym_matc8(matsymc8)
    if(sym.eq.'c2') call icos2f_sym_matc2(matsymc2)
    if(sym.eq.'icos')call icos2f_sym_mat(symmat)
    no0=-1
    do lz=-nx/2,nx/2-1
!**********************************************************************
        no1=int(float((lz)*100)/float(nx/2-1))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        do ky=-nx/2,nx/2-1
            do hx=-nx/2,nx/2-1
                r=sqrt(float(hx**2+ky**2+lz**2))
                if(r.gt.r0)  goto 99
                i3=(lz+nx/2)*nx**2+(ky+nx/2)*nx+(hx+nx/2)+1
                x0(1)=hx
                x0(2)=ky
                x0(3)=lz
                do isym=1,nsym
                  if(sym.eq.'icos')x=matmul(symmat(:,:,isym),x0)+nx/2
                  if(sym.eq.'c12')x=matmul(matsymc12(:,:,isym),x0)+nx/2
                  if(sym.eq.'c15')x=matmul(matsymc15(:,:,isym),x0)+nx/2
                  if(sym.eq.'c6')x=matmul(matsymc6(:,:,isym),x0)+nx/2
                  if(sym.eq.'c3')x=matmul(matsymc3(:,:,isym),x0)+nx/2
                  if(sym.eq.'c5a')x=matmul(matsymc5a(:,:,isym),x0)+nx/2
                  if(sym.eq.'c4')x=matmul(matsymc4(:,:,isym),x0)+nx/2
                  if(sym.eq.'c2')x=matmul(matsymc2(:,:,isym),x0)+nx/2
                  if(sym.eq.'c8')x=matmul(matsymc8(:,:,isym),x0)+nx/2
                   ix=int(x(1))
                   iy=int(x(2))
                   iz=int(x(3))
                   ax=x(1)-ix
                   ay=x(2)-iy
                   az=x(3)-iz
                   n1=iz*nx**2+iy*nx+ix+1
                   n2=n1+1
                   n3=n1+nx
                   n4=n3+1
                   n5=n1+nx**2
                   n6=n2+nx**2
                   n7=n3+nx**2
                   n8=n4+nx**2
                   img3d1(i3)=img3d1(i3)+ &
                             (1.0-ax)*(1.0-ay)*(1.0-az)*img3d0(n1)+ &
                             (    ax)*(1.0-ay)*(1.0-az)*img3d0(n2)+ &
                             (1.0-ax)*(    ay)*(1.0-az)*img3d0(n3)+ &
                             (    ax)*(    ay)*(1.0-az)*img3d0(n4)+ &
                             (1.0-ax)*(1.0-ay)*(    az)*img3d0(n5)+ &
                             (    ax)*(1.0-ay)*(    az)*img3d0(n6)+ &
                             (1.0-ax)*(    ay)*(    az)*img3d0(n7)+ &
                             (    ax)*(    ay)*(    az)*img3d0(n8)

                end do
99          end do
        end do
    end do

    do lz=-nx/2,nx/2-1
        do ky=-nx/2,nx/2-1
            do hx=-nx/2,nx/2-1
                r=sqrt(float(hx**2+ky**2+lz**2))
                if(r.gt.r0)  goto 991
                i3=(lz+nx/2)*nx**2+(ky+nx/2)*nx+(hx+nx/2)+1
!                if(lz.eq.0) then
!                    i3b=(lz+nx/2)*nx**2+(-ky+nx/2)*nx+(-hx+nx/2)+1
!                    img3d1(i3b)=img3d1(i3)
!                else
!                    ix=-hx;iy=-ky;iz=lz
!                    i3b=(iz+nx/2)*nx**2+(iy+nx/2)*nx+(ix+nx/2)+1
!                    img3d1(i3b)=img3d1(i3)
!                    i3c=(-lz+nx/2)*nx**2+(ky+nx/2)*nx+(-hx+nx/2)+1
!                    i3d=(-iz+nx/2)*nx**2+(iy+nx/2)*nx+(-ix+nx/2)+1
!                    img3d1(i3c)=img3d1(i3)
!                    img3d1(i3d)=img3d1(i3)
!                end if
991         end do
        end do
    end do

    img3d0=img3d1/float(nsym)
    deallocate(img3d1)
    return
end
