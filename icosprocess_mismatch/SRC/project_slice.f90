



subroutine project_slice(img3d,nx,theta0,phi0,omega0,slice)
    use common_icosproc
    implicit none
    integer nx,ix,iy,iz,i3,i2,n1,n2,n3,n4,n5,n6,n7,n8
    real theta0,phi0,omega0,theta,phi,omega
    real rotMat(3,3), Xprime(3), X(3), xbit, ybit, zbit
    real img3d(nx**3), slice(nx**2)

    open(10,file=trim(templateort),status='old')
    open(11,file=trim(model3d), form='unformatted',access='stream',status='old')
    open(21,file=trim(templatestck),form='unformatted',access='stream')

        theta=theta0
        phi  =phi0
        omega=omega0
        call euler_mat(theta,phi,omega,rotMat)

        do i3=1,nx**3
            iz=(i3-1)/nx**2
            iy=mod(i3-1,nx**2)/nx
            ix=mod(mod(i3-1,nx**2),nx)
            i2=iy*nx+ix+1
            Xprime(1)=float(ix)-nx/2
            Xprime(2)=float(iy)-nx/2
            Xprime(3)=float(iz)-nx/2
            if(sqrt(Xprime(1)**2+Xprime(2)**2+Xprime(3)**2).ge.(nx/2-1)) goto 99

            x=matmul(rotMat,Xprime)
            x=x+nx/2

            ix=int(X(1))
            iy=int(X(2))
            iz=int(X(3))

            xbit=X(1)-float(ix)
            ybit=X(2)-float(iy)
            zbit=X(3)-float(iz)

            n1= iz   *nx**2+ iy   *nx+ ix+1
            n2= iz   *nx**2+ iy   *nx+(ix+1)+1
            n3= iz   *nx**2+(iy+1)*nx+ ix+1
            n4= iz   *nx**2+(iy+1)*nx+(ix+1)+1
            n5=(iz+1)*nx**2+ iy   *nx+ ix+1
            n6=(iz+1)*nx**2+ iy   *nx+(ix+1)+1
            n7=(iz+1)*nx**2+(iy+1)*nx+ ix+1
            n8=(iz+1)*nx**2+(iy+1)*nx+(ix+1)+1


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

end subroutine

