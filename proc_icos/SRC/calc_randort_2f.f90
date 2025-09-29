subroutine calc_randort_2f
use common_proc_icos
implicit none
    integer particleID,imgID,imgparticleID,i,sumi
    real theta,phi,omega,centx0,centy0,PR,df1,df2,astig
    character(len=64):: char1
    real,allocatable::equ_ort(:,:),randx(:)
    allocate(equ_ort(3,60))

    open(11,file=trim(ort0),status='old')
    open(12,file=trim(newort))
    allocate(randx(100000))
    call randomreal(randx,100000)

    sumi=0
100 read(11,*,end=998) char1,particleID,theta,phi,omega,centx0,centy0,PR,imgID,imgparticleID,df1,df2,astig
        call equalort(theta,phi,omega,equ_ort,60)
!        do i=1,60
            sumi=sumi+1
            i=int(randx(sumi)*60.0)
!            print*,'sumi=',sumi,i
            write(12,900) char1,sumi,equ_ort(1,i),equ_ort(2,i),equ_ort(3,i),centx0,centy0,PR,imgID,imgparticleID,df1,df2,astig
!        end do
    goto 100
900 format(a64,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4)
998 print*,'the new ortfile is written in ',trim(newort)
999 end 
