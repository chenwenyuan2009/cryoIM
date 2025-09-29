subroutine equalort(theta,phi,omega,ort,N)
implicit none
    integer N,i
    real theta,phi,omega,rad2deg,ort(3,N),eularmat(3,3),r(3,3)
    real,allocatable::icos_rot(:,:,:)
    allocate(icos_rot(3,3,60))

    rad2deg=45.0/atan2(1.0,1.0)
    call Eularmatrix(theta,phi,omega,eularmat)
    call icos2f_operation_matrix(icos_rot,60)
    do i=1,60
        r=matmul(icos_rot(:,:,i),eularmat)
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

        

    