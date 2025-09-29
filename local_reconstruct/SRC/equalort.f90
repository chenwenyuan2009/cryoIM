subroutine equalort(theta,phi,omega,sym_rot3,ort,N)
implicit none
    integer N,i
    real theta,phi,omega,rad2deg,ort(3,N),eularmat(3,3),r(3,3),sym_rot3(3,3,N)
    character(len=64)::sym_matrix_file

    rad2deg=45.0/atan2(1.0,1.0)
    call Eularmatrix(theta,phi,omega,eularmat)
    do i=1,N
        r=matmul(sym_rot3(:,:,i),eularmat)
        if(r(3,3).ge.1.0) then
            ort(1,i)=0
        else if(r(3,3).lt.-1.0) then
            ort(1,i)=180.0
        else
            ort(1,i)=acos(r(3,3))*rad2deg
        end if
        ort(2,i)=atan2(r(2,3),r(1,3))*rad2deg
        ort(3,i)=atan2(r(3,2),-r(3,1))*rad2deg
    end do
    return
    end

        

    