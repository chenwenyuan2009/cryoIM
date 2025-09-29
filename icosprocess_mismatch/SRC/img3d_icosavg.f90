
subroutine img3d_icosavg(img3d0,nx,r0)
    implicit none
    integer nx,i3,nsym,hx,ky,lz,ncycle,ix,iy,iz,n1,n2,n3,n4,n5,n6,n7,n8,no0,no1
    integer i3b,i3c,i3d
    real img3d0(nx**3),symmat(3,3,60)
    real r,r0,x0(3),x(3),ax,ay,az
    real,allocatable::img3d1(:)
    allocate(img3d1(nx**3))

    print*,'performing real space icos average...'
    symmat=.0
    nsym  =1
    img3d1=.0
    nsym  =60
    call icos2f_sym_mat(symmat)

    no0=-1
    do lz=0,nx/2-1
!**********************************************************************
        no1=int(float((lz)*100)/float(nx/2-1))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        do ky=0,nx/2-1
            do hx=-nx/2,nx/2-1
                r=sqrt(float(hx**2+ky**2+lz**2))
                if(r.gt.r0)  goto 99
                i3=(lz+nx/2)*nx**2+(ky+nx/2)*nx+(hx+nx/2)+1
                x0(1)=hx
                x0(2)=ky
                x0(3)=lz
                do ncycle=1,nsym
                   x=matmul(symmat(:,:,ncycle),x0)+nx/2
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

    do lz=0,nx/2-1
        do ky=0,nx/2-1
            do hx=-nx/2,nx/2-1
                r=sqrt(float(hx**2+ky**2+lz**2))
                if(r.gt.r0)  goto 991
                i3=(lz+nx/2)*nx**2+(ky+nx/2)*nx+(hx+nx/2)+1
                if(lz.eq.0) then
                    i3b=(lz+nx/2)*nx**2+(-ky+nx/2)*nx+(-hx+nx/2)+1
                    img3d1(i3b)=img3d1(i3)
                else
                    ix=-hx;iy=-ky;iz=lz
                    i3b=(iz+nx/2)*nx**2+(iy+nx/2)*nx+(ix+nx/2)+1
                    img3d1(i3b)=img3d1(i3)
                    i3c=(-lz+nx/2)*nx**2+(ky+nx/2)*nx+(-hx+nx/2)+1
                    i3d=(-iz+nx/2)*nx**2+(iy+nx/2)*nx+(-ix+nx/2)+1
                    img3d1(i3c)=img3d1(i3)
                    img3d1(i3d)=img3d1(i3)
                end if
991         end do
        end do
    end do

    img3d0=img3d1/float(nsym)
    deallocate(img3d1)
    return
end
