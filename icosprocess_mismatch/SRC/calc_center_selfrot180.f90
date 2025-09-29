

subroutine calc_center_selfrot180
use common_icosproc
    implicit none
    integer pID,hx,ky,Nmax(1),ix,iy,ix0,iy0,i,j
    real a1(FFTsize*FFTsize),a2(FFTsize*FFTsize)
    real cc,centx,centy,z(-1:1,-1:1),dx,dy

    stckfile='abcdef'
    do pID=first,last
        if(trim(stckfile).ne.trim(particle(pID)%stckfile)) then
            close(11)
            stckfile=trim(particle(pID)%stckfile)
            open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(11) mrc
            call fseek(11,4*FFTsize**2*(particle(pID)%imgparticleID-1),1)
            read(11) a1
        else
            call fseek(11,4*FFTsize**2*(particle(pID)%imgparticleID-particle(pID-1)%imgparticleID-1),1)
            read(11) a1
        end if
        call img2d_lowpass(a1,FFTsize,maxR)
        if(proc2d%imgmask.eq.'y') call img2d_mask(a1,FFTsize,imgmask,0.0,0.0)
        call img2d_norm(a1,FFTsize)
        a2=.0
        !rotate 180 degree
        do i=1,FFTsize*FFTsize
            ky=(i-1)/FFTsize
            hx=mod(i-1,FFTsize)
            if((hx.eq.0).or.(ky.eq.0)) then
                a1(i)=.0
            else
                ky=FFTsize-ky
                hx=FFTsize-hx
                a2(ky*FFTsize+hx+1)=a1(i)
            end if
        end do
        call cross_correlation(a1,a2,FFTsize)

        Nmax=maxloc(a1)-1
        iy0=Nmax(1)/FFTsize
        ix0=mod(Nmax(1),FFTsize)
        do ky=-1,1
            do hx=-1,1
                ix=ix0+hx
                iy=iy0+ky
                j=iy*FFTsize+ix+1
                z(hx,ky)=a1(j)
            end do
        end do
        call guass_regression(z,dx,dy,cc)
        centx=(float(ix0-FFTsize/2)+dx)/2.0
        centy=(float(iy0-FFTsize/2)+dy)/2.0
        print*,pID,centx,centy,cc
        if(cc.ge. 1.0) cc=1.0
        if(cc.le.-1.0) cc=-1.0
        particle(pID)%x0 =particle(pID)%x
        particle(pID)%y0 =particle(pID)%y
        particle(pID)%PR0=particle(pID)%PR
        particle(pID)%x  =centx
        particle(pID)%y  =centy
        particle(pID)%PR =acos(cc)*45.0/atan2(1.0,1.0)
        particle(pID)%dx =particle(pID)%x -particle(pID)%x0
        particle(pID)%dy =particle(pID)%y -particle(pID)%y0
        particle(pID)%dPR=particle(pID)%PR-particle(pID)%PR0


    end do
    close(21)
    close(51)

    return
 end subroutine





