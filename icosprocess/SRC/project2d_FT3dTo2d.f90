
!根据取向，在三维傅立叶空间通过插值得到中央截面，然后通过傅立叶反变换得到二维投影.
subroutine project2d_FT3dTo2d(ft3d_r,ft3d_i,nx,theta0,phi0,omega0,FT2dr)
    implicit none
    integer nx,ix,iy,hx,ky,lz,n,i
    real rotMat(3,3),Xprime(3),XX(3),xbit,ybit,zbit
    real theta0,phi0,omega0,theta,phi,omega
    real tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
    real ft3d_r(nx**3),ft3d_i(nx**3),FT2dr(nx**2)
    real,allocatable::shift2d(:),FT2di(:)
    allocate(shift2d(nx**2),FT2di(nx**2))

    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        shift2d(i)=(-1)**(hx+ky)
    end do

    theta=theta0
    phi=phi0
    omega=omega0
    call euler_mat(theta,phi,omega,rotMat)

    n=0
    FT2dr=.0
    FT2di=.0
    do iy=-nx/2, nx/2-1
        do ix=-nx/2, nx/2-1
            n=n+1
            Xprime(1) = float(ix)
            Xprime(2) = float(iy)
            Xprime(3) = .0
            xx=matmul((rotMat),Xprime)

            if( sqrt(XX(1)**2+XX(2)**2+XX(3)**2) .gt. nx/2-2) goto 99
            XX(1) = XX(1) + nx/2
            XX(2) = XX(2) + nx/2
            XX(3) = XX(3) + nx/2

            hx = int(XX(1))
            ky = int(XX(2))
            lz = int(XX(3))
            xbit = XX(1) - float(hx)
            ybit = XX(2) - float(ky)
            zbit = XX(3) - float(lz)

            tmp1=ft3d_r(lz    *nx**2+ ky   *nx+hx+1)
            tmp2=ft3d_r(lz    *nx**2+ ky   *nx+hx+1 +1)
            tmp3=ft3d_r(lz    *nx**2+(ky+1)*nx+hx+1)
            tmp4=ft3d_r(lz    *nx**2+(ky+1)*nx+hx+1 +1)
            tmp5=ft3d_r((lz+1)*nx**2+ ky   *nx+hx+1)
            tmp6=ft3d_r((lz+1)*nx**2+ ky   *nx+hx+1 +1)
            tmp7=ft3d_r((lz+1)*nx**2+(ky+1)*nx+hx+1)
            tmp8=ft3d_r((lz+1)*nx**2+(ky+1)*nx+hx+1 +1)

            FT2dr(n) =   &
                         (1.-xbit)*(1.-ybit)*(1.-zbit) * tmp1 +  &
                             xbit *(1.-ybit)*(1.-zbit) * tmp2 +  &
                         (1.-xbit)*    ybit *(1.-zbit) * tmp3 +  &
                             xbit *    ybit *(1.-zbit) * tmp4 +  &
                         (1.-xbit)*(1.-ybit)*    zbit  * tmp5 +  &
                             xbit *(1.-ybit)*    zbit  * tmp6 +  &
                         (1.-xbit)*    ybit *    zbit  * tmp7 +  &
                             xbit *    ybit *    zbit  * tmp8

                tmp1=ft3d_i(lz    *nx**2+ ky   *nx+hx+1)
                tmp2=ft3d_i(lz    *nx**2+ ky   *nx+hx+1 +1)
                tmp3=ft3d_i(lz    *nx**2+(ky+1)*nx+hx+1)
                tmp4=ft3d_i(lz    *nx**2+(ky+1)*nx+hx+1 +1)
                tmp5=ft3d_i((lz+1)*nx**2+ ky   *nx+hx+1)
                tmp6=ft3d_i((lz+1)*nx**2+ ky   *nx+hx+1 +1)
                tmp7=ft3d_i((lz+1)*nx**2+(ky+1)*nx+hx+1)
                tmp8=ft3d_i((lz+1)*nx**2+(ky+1)*nx+hx+1 +1)

            FT2di(n) =   &
                         (1.-xbit)*(1.-ybit)*(1.-zbit) * tmp1 +  &
                             xbit *(1.-ybit)*(1.-zbit) * tmp2 +  &
                         (1.-xbit)*    ybit *(1.-zbit) * tmp3 +  &
                             xbit *    ybit *(1.-zbit) * tmp4 +  &
                         (1.-xbit)*(1.-ybit)*    zbit  * tmp5 +  &
                             xbit *(1.-ybit)*    zbit  * tmp6 +  &
                         (1.-xbit)*    ybit *    zbit  * tmp7 +  &
                             xbit *    ybit *    zbit  * tmp8
99      end do
    end do
    FT2dr=FT2dr*shift2d
    FT2di=FT2di*shift2d
    call img2d_FFT(FT2dr,FT2di,nx,nx,-1)
    FT2dr=FT2dr*shift2d

    deallocate(shift2d,FT2di)
end subroutine


