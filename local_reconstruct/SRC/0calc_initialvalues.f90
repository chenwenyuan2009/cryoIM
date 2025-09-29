subroutine calc_initialvalues
use common_local_reconstruct
    implicit none
    integer particleID,imgID,imgparticleID
    real theta,phi,omega,x,y,PR,df1,df2,astigang
    character(len=64)::char

    if(proc2d%imgstck.eq.'y') then
        open(1,file=trim(imgstck),form='unformatted',access='stream',status='old')
        read(1) mrc
        close(1)
    else
        open(1,file=trim(ortfile),status='old')
        read(1,*) char
        close(1)
        open(2,file=trim(char),form='unformatted',access='stream',status='old')
        read(2) mrc
        close(2)
    end if
    FFTsize0=mrc%nx
    FFTsize=FFTsize
    if(proc2d%imgmask.eq.'y') then
        imgmask=imgmask/apix
        if(imgmask.ge.float(FFTsize)/2.0) imgmask=float(FFTsize)/2.0-1
    else
        imgmask=float(FFTsize)/2.0-1
    end if
    maxRes=apix*FFTsize/maxRes
    maxR  =int(maxRes)
    if(maxR.ge.FFTsize/2) maxR=FFTsize/2-2
    !if(mod(maxR,2).ne.0) maxR=maxR+1

    if((trim(sym).eq.'I2').or.(trim(sym).eq.'i2')) Nsym=60
    if((trim(sym).eq.'I3').or.(trim(sym).eq.'i3')) Nsym=60
    if((trim(sym).eq.'I5').or.(trim(sym).eq.'i5')) Nsym=60
    if((trim(sym).eq.'c1').or.(trim(sym).eq.'C1')) Nsym=1
    if((trim(sym).eq.'c2').or.(trim(sym).eq.'C2')) Nsym=2
    if((trim(sym).eq.'c3').or.(trim(sym).eq.'C3')) Nsym=3
    if((trim(sym).eq.'c4').or.(trim(sym).eq.'C4')) Nsym=4
    if((trim(sym).eq.'c5').or.(trim(sym).eq.'C5')) Nsym=5
    if((trim(sym).eq.'c6').or.(trim(sym).eq.'C6')) Nsym=6
    if((trim(sym).eq.'c7').or.(trim(sym).eq.'C7')) Nsym=7
    if((trim(sym).eq.'c8').or.(trim(sym).eq.'C8')) Nsym=8
    if((trim(sym).eq.'c9').or.(trim(sym).eq.'C9')) Nsym=9
    if((trim(sym).eq.'c10').or.(trim(sym).eq.'C_10')) Nsym=10
    if((trim(sym).eq.'c11').or.(trim(sym).eq.'C_11')) Nsym=11
    if((trim(sym).eq.'c12').or.(trim(sym).eq.'C_12')) Nsym=12
    if((trim(sym).eq.'c13').or.(trim(sym).eq.'C_13')) Nsym=13
    if((trim(sym).eq.'c14').or.(trim(sym).eq.'C_14')) Nsym=14
    if((trim(sym).eq.'c15').or.(trim(sym).eq.'C_15')) Nsym=15
    allocate(matsym(3,3,Nsym))
    matsym=.0
    matsym(1,1,:)=1;matsym(2,2,:)=1;matsym(3,3,:)=1
    if((trim(sym).eq.'i2').or.(trim(sym).eq.'I2')) then
        call icos_matrix_i2(matsym,Nsym)
    else if((trim(sym).eq.'i3').or.(trim(sym).eq.'I3')) then
        call icos_matrix_i3(matsym,Nsym)
    else if((trim(sym).eq.'i5').or.(trim(sym).eq.'I5')) then
        call icos_matrix_i5(matsym,Nsym)
    else
        call C_matrix(matsym,Nsym)
    end if

    totalparticles=0
    open(11,file=trim(ortfile),status='old')
100 read(11,*,end=900) char,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
        totalparticles=totalparticles+1
        goto 100
900 close(11)
    if(proc2d%last.eq.'n') last=totalparticles
    if((proc2d%last.eq.'y').and.(last.gt.totalparticles)) last=totalparticles
    allocate(particle(first:last))
    particle(:)%quality='b'
    allocate(ctf2d(FFTsize**2))
    ctf2d=1.0

    return
end

