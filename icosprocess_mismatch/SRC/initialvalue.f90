!各参数赋初值
subroutine initialvalue
use common_icosproc
    implicit none
    character(len=128)::tmp

    if(proc2d%imgstck.eq.'y') then
        open(2,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(2) mrc
        close(2)
    else
       if (proc2d%ortfile.eq.'y')then
        open(1,file=trim(ortfile),status='old')
        read(1,*) tmp
        close(1)
        open(2,file=trim(tmp),form='unformatted',access='stream',status='old')
        read(2) mrc
        close(2)
        else
        open(2,file=trim(model3d),form='unformatted',access='stream',status='old')
        read(2) mrc
        close(2)
        endif
    end if

    if(proc2d%ortfile.eq.'y')then
    open(1,file=trim(ortfile),status='old')
    last=0
100 read(1,*,end=110) tmp
    last=last+1
    goto 100
110 close(1)
    endif
    FFTsize       =mrc%nx
    geneFFTsize =FFTsize
    inter_box     =2
    inter_w       =1
    PR_threshold  =-1.0
    centbound     =3
    CTFpar%vol    =300
    CTFpar%Cs     =2.7
    CTFpar%ampwg  =0.07
    CTFpar%Bfactor=.0
    first         =1
    minRes        =300
    maxRes        =8.0
    apix          =1.0
    Nsym          =1
    Nlin          =1
    searchstep    =2.0
    Ncycle        =0
    Padsize       =1
    result3d='map3d.mrc'

    return
end
