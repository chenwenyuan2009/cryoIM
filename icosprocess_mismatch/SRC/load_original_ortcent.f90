subroutine load_original_ortcent
use common_icosproc
    implicit none
    integer i
    open(11,file=trim(ortfile),status='old')

    particle(:)%dtheta=.0;particle(:)%dphi=.0;particle(:)%domega=.0;particle(:)%dx=.0;particle(:)%dy=.0;particle(:)%dPR=.0
    useparticle=0
    do i=1,first-1
        read(11,*)
    end do
    do i=first,last
        read(11,*) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                   particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                   particle(i)%df1,particle(i)%df2,particle(i)%astigang
        if((particle(i)%PR.ge.PR_threshold).and.(sqrt(particle(i)%x**2+particle(i)%y**2).le.centbound)&
                   .and.(particle(i)%PR.lt.PR_threshold_h)) then
                   useparticle=useparticle+1
        endif
    end do
    close(11)
    return
end subroutine


