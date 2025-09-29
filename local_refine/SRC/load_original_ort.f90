subroutine load_original_ort
use common_loc_refine
    implicit none
    integer i
    open(11,file=trim(local_ort0),status='old')

    do i=1,first-1
       read(11,*)
    end do
    do i=first,last
  
        read(11,*)particle(i)%stckfile, particle(i)%particleID, particle(i)%theta, particle(i)%phi, particle(i)%omega, &
                  particle(i)%x, particle(i)%y, particle(i)%PR, particle(i)%imgID, particle(i)%imgparticleID,  &
                  particle(i)%df1, particle(i)%df2, particle(i)%astigang, &
                  particle(i)%local_imgID, particle(i)%local_hx0, particle(i)%local_ky0, particle(i)%local_theta0, particle(i)%local_phi0, &
                  particle(i)%local_omega0, particle(i)%local_x0, particle(i)%local_y0, particle(i)%local_z0,particle(i)%local_PR0
        particle(i)%quality='g'
    end do
    close(11)
    return
end subroutine


