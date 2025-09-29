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
    if(mode.ne.6)then
        read(11,*) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                   particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                   particle(i)%df1,particle(i)%df2,particle(i)%astigang
    else
         read(11,*) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                   particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                   particle(i)%df1,particle(i)%df2,particle(i)%astigang,particle(i)%x0,particle(i)%y0
                   
    endif
        if(particle(i)%PR.le.PR_threshold.and.(abs(particle(i)%x).le.boundX ).and.(abs(particle(i)%y).le.boundX))then
         useparticle=useparticle+1
         endif
    end do
    close(11)
    return
end subroutine


