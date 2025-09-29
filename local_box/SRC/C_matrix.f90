subroutine C_matrix(r,Nsym)
    implicit none
    integer i1,Nsym
    real pi,phi
    real r(3,3,Nsym)

    pi=4.0*atan2(1.0,1.0)
    phi=2.0*pi/float(Nsym)
    r(1,1,1)=cos(phi); r(1,2,1)=-sin(phi); r(1,3,1)= .0
    r(2,1,1)=sin(phi); r(2,2,1)= cos(phi); r(2,3,1)= .0
    r(3,1,1)=.0      ; r(3,2,1)=.0       ; r(3,3,1)=1.0

    do i1=2,Nsym
        r(:,:,i1)=matmul(r(:,:,i1-1),r(:,:,1))
    end do
    return
end
