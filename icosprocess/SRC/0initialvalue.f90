! initial parameters
subroutine initialvalue
use common_icosproc
    implicit none
    character(len=128)::tmp

    if(proc2d%imgstck.eq.'y') then
        open(1,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(1) mrc
        close(1)
        open(1,file=trim(ortfile),status='old')
    else
        open(1,file=trim(ortfile),status='old')
        read(1,*) tmp
        rewind(1)
        open(2,file=trim(tmp),form='unformatted',access='stream',status='old')
        read(2) mrc
        close(2)
    end if

    last=0
100 read(1,*,end=110) tmp
    last=last+1
    goto 100
110 close(1)
    FFTsize       =mrc%nx
    imgmask=FFTsize/2-2
    boundX= FFTsize/2-2
    inter_box     =2
    inter_w       =1
    PR_threshold  =90.0
    CTFpar%vol    =300
    CTFpar%Cs     =2.7
    CTFpar%ampwg  =0.1
    CTFpar%Bfactor=.0
    first         =1
    minRes        =300
    maxRes        =8.0
    apix          =1.0
    Nsym          =60
    Nlin          =60
    Ntemplate     =20
    searchstep    =2.0
    Ncycle        =6
    loc_box       =FFTsize
    result3d='map3d.mrc'

    return
end
