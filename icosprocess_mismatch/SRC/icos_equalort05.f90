subroutine icos_equalort05(theta0,phi0,omega0,angle)
    implicit none
    integer i
    real theta0,phi0,omega0,theta,phi,omega,angle(3,5),maticos(3,3,5),mat(3,3),matoular(3,3)

    theta=theta0
    phi  =-phi0
    omega=180.0-omega0
    call icos2f_sym_mat05(maticos)
    call euler_mat(theta,phi,omega,matoular)

    do i=1,5
        mat=matmul(maticos(:,:,i),matoular)
        angle(1,i)=acos(mat(3,3))*45.0/atan2(1.0,1.0)
        angle(2,i)=-atan2(mat(2,3),mat(1,3))*45.0/atan2(1.0,1.0)
        angle(3,i)=atan2(mat(3,2),mat(3,1))*45.0/atan2(1.0,1.0)
    end do
    return
end
