subroutine calc_initialvilues
use common_loc_refine
    implicit none

    integer particleID,imgID,imgparticleID
    real theta,phi,omega,x,y,PR,df1,df2,astigang
    character(len=64)::char

    open(11,file=trim(model3d),form='unformatted',access='stream',status='old')
    read(11) mrc
    close(11) 
    FFTsize=mrc%nx

    if((trim(sym).eq.'c1').or.(trim(sym).eq.'C1')) Nsym=1
    if((trim(sym).eq.'c2').or.(trim(sym).eq.'C2')) Nsym=2
    if((trim(sym).eq.'c3').or.(trim(sym).eq.'C3')) Nsym=3
    if((trim(sym).eq.'c4').or.(trim(sym).eq.'C4')) Nsym=4
    if((trim(sym).eq.'c5').or.(trim(sym).eq.'C5')) Nsym=5
    if((trim(sym).eq.'c6').or.(trim(sym).eq.'C6')) Nsym=6
    if((trim(sym).eq.'c7').or.(trim(sym).eq.'C7')) Nsym=7
    if((trim(sym).eq.'c8').or.(trim(sym).eq.'C8')) Nsym=8
    if((trim(sym).eq.'c9').or.(trim(sym).eq.'C9')) Nsym=9
    if((trim(sym).eq.'c10').or.(trim(sym).eq.'C10')) Nsym=10
    if((trim(sym).eq.'c11').or.(trim(sym).eq.'C11')) Nsym=11
    if((trim(sym).eq.'c12').or.(trim(sym).eq.'C12')) Nsym=12
    if((trim(sym).eq.'c13').or.(trim(sym).eq.'C13')) Nsym=13
    if((trim(sym).eq.'c14').or.(trim(sym).eq.'C14')) Nsym=14
    if((trim(sym).eq.'c15').or.(trim(sym).eq.'C15')) Nsym=15
    if((trim(sym).eq.'asym').or.(trim(sym).eq.'ASYM')) Nsym=180
    if((trim(sym).eq.'c180').or.(trim(sym).eq.'C180')) Nsym=180
    if((trim(sym).eq.'c360').or.(trim(sym).eq.'C360')) Nsym=360
    if((trim(sym).eq.'i2').or.(trim(sym).eq.'I2')) Nsym=60
    if((trim(sym).eq.'i5').or.(trim(sym).eq.'I5')) Nsym=60

    if(proc2d%imgstck.eq.'y') then
        open(21,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(21) mrc
        close(21)
    else
        open(21,file=trim(local_ort0),status='old')
        read(21,*) char
        open(31,file=trim(char),form='unformatted',access='stream',status='old')
        read(31) mrc
        close(21)
        close(31)
    end if
    FFTsize0=mrc%nx

    totalparticles=0
    open(21,file=trim(local_ort0),status='old')
100 read(21,*,end=900) char,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
        totalparticles=totalparticles+1
        goto 100
900 close(21)
    if(proc2d%last.eq.'n') last=totalparticles
    if((proc2d%last.eq.'y').and.(last.gt.totalparticles)) last=totalparticles
    allocate(particle(first:last))
    particle(:)%quality='b'
    call load_original_ort

    if(proc2d%imgmask.eq.'y')  then
        imgmask=imgmask/apix
    else
        imgmask=FFTsize/2-2
    end if

    maxRes=apix*real(FFTsize)/maxRes
    minRes=apix*real(FFTsize)/minRes

    if(proc2d%refine_angle.eq.'n') refine_angle=360.0/(maxRes*sqrt(5.0)*3.14159); !angle step for ort refine
    if(proc2d%refine_shift.eq.'n') refine_shift=float(FFTsize)/(maxRes*sqrt(5.0)); !cent shift step for cent refine
    return
end

