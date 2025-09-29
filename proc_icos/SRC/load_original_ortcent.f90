
subroutine load_original_ortcent
use common_proc_icos
    implicit none
    integer particleID,imgID,imgparticleID,i
    real theta,phi,omega,centx0,centy0,PR,df1,df2,astig
    character(len=64)::char1
    open(11,file=trim(ort0),status='old')

    totalparticle=0
100 read(11,*,end=200) char1,particleID,theta,phi,omega,centx0,centy0,PR,imgID,imgparticleID,df1,df2,astig
    totalparticle=totalparticle+1
    goto 100
200 rewind(11)
    allocate(particle(totalparticle))
    
    particle(:)%dtheta=.0;particle(:)%dphi=.0;particle(:)%domega=.0;particle(:)%dx=.0;particle(:)%dy=.0;particle(:)%dPR=.0
    particle(:)%status='g'
    do i=1,totalparticle
            read(11,*) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                       particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                       particle(i)%df1,particle(i)%df2,particle(i)%astigang
            if(particle(i)%PR.gt.PR_threshold) particle(i)%status='b'
    end do
    close(11)
    return
end subroutine


