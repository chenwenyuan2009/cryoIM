subroutine checkstckfile
use common_icosproc
    implicit none
    integer i,n

    n=0
    do i=first,last
        open(21,file=trim(particle(i)%stckfile),form='unformatted',access='stream',status='old',err=900)
        read(21) mrc
        if(mrc%nz.lt.particle(i)%imgparticleID) then
            print*,'the particle number between ortfile and ',trim(particle(i)%stckfile),' is not match.'
            n=n+1
        end if
        close(21)
        goto 910
900     print*,'the stck file ',trim(particle(i)%stckfile),' is not exist.'
910 end do

    if(n.ne.0) then
        print*,'please check the ortfile and stckfiles.'
        stop
    else
        print*,'the ortfile and stckfiles are match each other.'
    end if
    return
end


