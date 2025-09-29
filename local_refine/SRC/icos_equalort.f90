
subroutine icos_equalort(theta0,phi0,omega0,angle)
   use common_loc_refine
    implicit none
    integer i
    real theta0,phi0,omega0,theta,phi,omega,angle(3,Nsym),maticos(3,3,Nsym),mat(3,3),matoular(3,3)

    theta=theta0
    phi  =-phi0
    omega=180.0-omega0
    if(sym.eq.'i2'.or.sym.eq.'I2')then 
    call icos_matrix_i2(maticos,Nsym)
    elseif(sym.eq.'i5'.or.sym.eq.'I5')then
    call icos_matrix_i5(maticos,Nsym)
    else
    call icos2f_sym_matc(maticos,Nsym)
    endif
    call euler_mat(theta,phi,omega,matoular)

    do i=1,Nsym
        mat=matmul(maticos(:,:,i),matoular)
        angle(1,i)=acos(mat(3,3))*45.0/atan2(1.0,1.0)
        angle(2,i)=-atan2(mat(2,3),mat(1,3))*45.0/atan2(1.0,1.0)
        angle(3,i)=atan2(mat(3,2),mat(3,1))*45.0/atan2(1.0,1.0)
    end do
    return
end

subroutine icos_equalortasym(theta0,phi0,omega0,angle)
   use common_loc_refine
    implicit none
    integer i
    real theta0,phi0,omega0,theta,phi,omega,angle(3,Nsym),maticos(3,3,Nsym),mat(3,3),matoular(3,3)

    theta=theta0
    phi  =-phi0
    omega=180.0-omega0
    call icos2f_sym_matasym(maticos,Nsym)
    call euler_mat(theta,phi,omega,matoular)

    do i=1,Nsym
        mat=matmul(maticos(:,:,i),matoular)
        angle(1,i)=acos(mat(3,3))*45.0/atan2(1.0,1.0)
        angle(2,i)=-atan2(mat(2,3),mat(1,3))*45.0/atan2(1.0,1.0)
        angle(3,i)=atan2(mat(3,2),mat(3,1))*45.0/atan2(1.0,1.0)
    end do
    return
end
