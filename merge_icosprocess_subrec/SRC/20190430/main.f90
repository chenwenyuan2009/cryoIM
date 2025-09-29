
program main
    use common_merge_icosprocess_subrec
    implicit none
    integer iargc,n,i


    n=iargc()
    if (n.lt.2) then
      call help
      stop
    end if
    do i=1,n
        call getarg(i,par(i))
    end do
    nfile=n
    icosavg='n'
    open(11,file=trim(par(1))//'oddtnf',form='unformatted',access='stream',status='old',err=900)
    goto 901
900 open(11,file=trim(par(1)),form='unformatted',access='stream',status='old')
901 read(11) FFTsize,maxR,apix,imgmask,nFSC,nicosavg
    print*,FFTsize,maxR,apix,imgmask,nFSC,nicosavg
    if(nicosavg.eq..1) icosavg='y'
    close(11)
    if(nFSC.eq.0) then
        call merge_reconstruct
    else if(nFSC.eq.1) then
        call merge_reconstruct_FSC
    end if

    end


