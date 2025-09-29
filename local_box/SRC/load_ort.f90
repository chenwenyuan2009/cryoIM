subroutine load_ort
use common_icos_loc_box
    implicit none
    integer i,particleID,imgID,imgparticleID
    real theta,phi,omega,x0,y0,PR,df1,df2,dfangle
    character(len=64)::stckfile
    open(11,file=trim(icos_ort),status='old')

    totalparticle=0
100 read(11,*,end=110) stckfile,particleID,theta,phi,omega,x0,y0,PR,imgID,imgparticleID,df1,df1,dfangle
    totalparticle=totalparticle+1
 
    goto 100
110 rewind(11)
    allocate(particle(totalparticle))
    print*,'totalparticle=',totalparticle
    if(proc2d_last.eq.'n') last=totalparticle
    if(proc2d_last.eq.'y') then
        if(last.gt.totalparticle) then
            print*, 'last should smaller than totalparticle'
            stop
        end if
    end if

    do i=1,totalparticle
        read(11,*) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                   particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                   particle(i)%df1,particle(i)%df2,particle(i)%astigang
        if(particle(i)%PR.le.PR_threshold .and. abs(particle(i)%x).le.boundX .and. abs(particle(i)%y).le.boundX) then
            particle(i)%status='g'
        else
            particle(i)%status='b'
        end if
    end do
    close(11)

    open(12,file=trim(particle(1)%stckfile),form='unformatted',access='stream',status='old')
    read(12) mrc
    FFTsize=mrc%nx
    close(12)
    return
end subroutine


