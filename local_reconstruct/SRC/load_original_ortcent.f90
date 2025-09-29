subroutine load_original_ortcent
use common_local_reconstruct
    implicit none
    integer i
    open(11,file=trim(ortfile),status='old')

    useparticle=0
    do i=1,first-1
       read(11,*)
    end do
    do i=first,last
        read(11,*)particle(i)%stckfile, particle(i)%particleID, particle(i)%theta, particle(i)%phi, particle(i)%omega, &
                  particle(i)%x, particle(i)%y, particle(i)%PR, particle(i)%imgID, particle(i)%imgparticleID,  &
                  particle(i)%df1, particle(i)%df2, particle(i)%astigang, &
                  particle(i)%local_imgID, particle(i)%local_hx0, particle(i)%local_ky0, particle(i)%local_theta, particle(i)%local_phi, &
                  particle(i)%local_omega, particle(i)%local_x, particle(i)%local_y, particle(i)%local_z

        if(particle(i)%PR.le.PR_threshold.and.abs(particle(i)%local_x).lt.boundX.and.abs(particle(i)%local_y).lt.boundX)then
            if(abs(particle(i)%local_theta-90.0).lt.tilt)then
            useparticle=useparticle+1
            particle(i)%quality='g'
            endif
        end if
    end do
    close(11)
    return
end subroutine


