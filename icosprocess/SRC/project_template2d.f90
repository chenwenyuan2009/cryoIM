

subroutine project_template2d
    use common_icosproc
    implicit none
    integer ix,iy,iz,i3,i2,itemplate,n1,n2,n3,n4,n5,n6,n7,n8
    real theta,phi,omega
    real rotMat(3,3), Xprime(3), X(3), xbit, ybit, zbit
    real,allocatable :: img3d(:), slice(:)

    allocate(img3d(FFTsize**3))
    allocate(slice(FFTsize**2))

    open(11,file=trim(model3d), form='unformatted',access='stream',status='old')
!    open(21,file='template.mrc',form='unformatted',access='stream')

    img3d     =.0
    template2d=.0
    read(11) mrc
    read(11) img3d
    close(11)
!    mrc%nz=totaltemplate
!    write(21) mrc

    do itemplate=1, TotalTemplate
        slice=.0
        theta=template(itemplate)%theta
        phi  =template(itemplate)%phi
        omega=template(itemplate)%omega
        print*,itemplate,theta,phi,omega
        call euler_mat(theta,phi,omega,rotMat)

        do i3=1,FFTsize**3
            iz=(i3-1)/FFTsize**2
            iy=mod(i3-1,FFTsize**2)/FFTsize
            ix=mod(mod(i3-1,FFTsize**2),FFTsize)
            i2=iy*FFTsize+ix+1
            Xprime(1)=float(ix)-FFTsize/2
            Xprime(2)=float(iy)-FFTsize/2
            Xprime(3)=float(iz)-FFTsize/2
            if(sqrt(Xprime(1)**2+Xprime(2)**2+Xprime(3)**2).ge.(FFTsize/2-1)) goto 99

            x=matmul(rotMat,Xprime)
            x=x+FFTsize/2

            ix=int(X(1))
            iy=int(X(2))
            iz=int(X(3))

            xbit=X(1)-float(ix)
            ybit=X(2)-float(iy)
            zbit=X(3)-float(iz)

            n1= iz   *FFTsize**2+ iy   *FFTsize+ ix+1
            n2= iz   *FFTsize**2+ iy   *FFTsize+(ix+1)+1
            n3= iz   *FFTsize**2+(iy+1)*FFTsize+ ix+1
            n4= iz   *FFTsize**2+(iy+1)*FFTsize+(ix+1)+1
            n5=(iz+1)*FFTsize**2+ iy   *FFTsize+ ix+1
            n6=(iz+1)*FFTsize**2+ iy   *FFTsize+(ix+1)+1
            n7=(iz+1)*FFTsize**2+(iy+1)*FFTsize+ ix+1
            n8=(iz+1)*FFTsize**2+(iy+1)*FFTsize+(ix+1)+1


            slice(i2) = slice(i2)+  &
                                  (1.-xbit)*(1.-ybit)*(1.-zbit) * img3d(n1) +  &
                                      xbit *(1.-ybit)*(1.-zbit) * img3d(n2) +  &
                                  (1.-xbit)*    ybit *(1.-zbit) * img3d(n3) +  &
                                      xbit *    ybit *(1.-zbit) * img3d(n4) +  &
                                  (1.-xbit)*(1.-ybit)*    zbit  * img3d(n5) +  &
                                      xbit *(1.-ybit)*    zbit  * img3d(n6) +  &
                                  (1.-xbit)*    ybit *    zbit  * img3d(n7) +  &
                                      xbit *    ybit *    zbit  * img3d(n8)
99      end do
!        write(21) slice
        template2d(:,1,itemplate)=slice(:)
    end do
    deallocate(img3d)
    deallocate(slice)

end subroutine

