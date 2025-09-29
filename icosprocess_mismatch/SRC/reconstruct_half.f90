


subroutine reconstruct_half
use common_icosproc
    implicit none
    integer i,i2,j,hx,ky,imgID0,no0,no1,imgparticle0,tmp1
    integer usedparticle,signz,tnfsize
    integer isym,lx(3),cubID(8)
    real theta,phi,omega,dr,rot(3,3),x0(3),x(3),x00(3),boundX
    real matsym(3,3,60),ax(3),disx(8)
    real,allocatable::a(:),b(:),shift2d(:)
    real,allocatable::fa(:),fb(:),wg(:),img2dstck(:)

    tnfsize=(maxR+1)*(2*maxR+1)**2
    allocate(a(FFTsize**2),b(FFTsize**2),shift2d(FFTsize**2))
    allocate(fa(tnfsize))
    allocate(fb(tnfsize))
    allocate(wg(tnfsize))
    boundX=FFTsize/4.0
    fa=.0;fb=.0;wg=.0
    imgID0=-100
    imgparticle0=-100
    usedparticle=0
    matsym=.0
    matsym(1,1,1)=1;matsym(2,2,1)=1;matsym(3,3,1)=1
    if(sym.eq.'icos') call icos2f_sym_mat(matsym)

    do i2=1,FFTsize**2
        ky=(i2-1)/FFTsize
        hx=mod(i2-1,FFTsize)
        shift2d(i2)=(-1)**(hx+ky)
    end do

    no0=-1
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2))
    print*,'perform half particle...'
    do i=first,last,2
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first-2))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(11,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(11) mrc
            allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))
            read(11) img2dstck
            close(11)
        end if
        a(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        if((particle(i)%PR.le.PR_threshold).and.(abs(particle(i)%x).le.boundX ).and.(abs(particle(i)%y).le.boundX)) then
            usedparticle=usedparticle+1
            call img2d_norm(a,FFTsize)
            if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
                if(particle(i)%imgID.ne.imgID0) then
                    call getCTF2d(FFTsize,particle(i)%df1,particle(i)%df2,particle(i)%astigang,  &
                                apix,CTFpar%vol,CTFpar%Cs,CTFpar%ampwg,ctf2d)
                    imgID0=particle(i)%imgID
                end if
            end if

            if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a,FFTsize,imgmask,particle%x,particle%y)
            a=a*shift2d
            b=.0
            call img2d_FFT(a,b,FFTsize,FFTsize,1)
            call img2d_FT_centshift(a,b,FFTsize,particle(i)%x,particle(i)%y)

            theta=particle(i)%theta
            phi  =particle(i)%phi
            omega=particle(i)%omega
            call euler_mat(theta,phi,omega,rot)

            do ky=-maxR,maxR
                do hx=0,maxR
                    dr=sqrt(float(hx**2+ky**2))
                    if(dr.le.maxR) then
                        i2=(ky+FFTsize/2)*FFTsize+(hx+FFTsize/2)+1
                        x(1)=hx
                        x(2)=ky
                        x(3)=.0
                        x0=matmul(rot,x)
!************************************************************************************
                        x00=x0
                        do isym=1,Nsym
                            x=x00
                            x0=matmul(matsym(:,:,isym),x)
                            if(x0(3).lt..0) then
                                x0=-x0
                                signz=-1
                            else
                                signz=1
                            end if
                            x0(1:2)=x0(1:2)+maxR
!************************************************************************************
                            lx=int(x0)
                            ax=x0-lx
                            disx(1)=(1.0-ax(1))*(1.0-ax(2))*(1.0-ax(3))
                            disx(2)=(    ax(1))*(1.0-ax(2))*(1.0-ax(3))
                            disx(3)=(1.0-ax(1))*(    ax(2))*(1.0-ax(3))
                            disx(4)=(    ax(1))*(    ax(2))*(1.0-ax(3))
                            disx(5)=(1.0-ax(1))*(1.0-ax(2))*(    ax(3))
                            disx(6)=(    ax(1))*(1.0-ax(2))*(    ax(3))
                            disx(7)=(1.0-ax(1))*(    ax(2))*(    ax(3))
                            disx(8)=(    ax(1))*(    ax(2))*(    ax(3))
                            cubID(1)=lx(3)*(2*maxR+1)**2+lx(2)*(2*maxR+1)+lx(1)+1
                            cubID(2)=cubID(1)+1
                            cubID(3:4)=cubID(1:2)+(2*maxR+1)
                            cubID(5:8)=cubID(1:4)+(2*maxR+1)**2
                            do j=1,8
                                fa(cubID(j))=fa(cubID(j))+a(i2)*disx(j)*ctf2d(i2)
                                fb(cubID(j))=fb(cubID(j))+b(i2)*disx(j)*ctf2d(i2)*signz
                                wg(cubID(j))=wg(cubID(j))+      disx(j)*ctf2d(i2)**2
                            end do
                        end do
                    end if
                end do
            end do
        end if
    end do
    print*,'writing tnf data...'
    open(21,file=trim(tnffile), form='unformatted',access='stream',status='replace')
    tmp1=0
    if((proc3d%subrec.eq.'y').or.(proc3d%subrec.eq.'Y')) tmp1=1
    write(21) FFTsize,maxR,apix,imgmask,tmp1
    write(21) fa
    write(21) fb
    write(21) wg
    close(21)
    deallocate(fa,fb,wg)
    return
end



