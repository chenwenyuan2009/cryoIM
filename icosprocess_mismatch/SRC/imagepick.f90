subroutine imagepick
    use common_icosproc
     implicit none
    integer i,j,imgID0,no0,no1
    real Pr
    real,allocatable::a0(:),img2dstck(:)

    allocate(a0(FFTsize**2),img2dstck(FFTsize**2))

    open(21,file=trim(newortfile))

    no0=-1
    imgID0=-100
    !read(12) mrc
    stckfile='abcde'
    do i=first,last
!**********************************************************************
        no1=int(float((i-first)*100)/float(last-first))
        if((no1.gt.no0).and.(mod(no1,10).eq.0)) then
            no0=no1
            write(*,101) no0,'%'
        end if
101     FORMAT(i7,a1,$)
!**********************************************************************
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(12,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(12) mrc
            allocate(img2dstck(FFTsize**2*mrc%nz))
            read(12)img2dstck
            close(12)
        end if
        a0(:)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+1:particle(i)%imgparticleID*FFTsize**2)
        call img2d_centshift_FT(a0,FFTsize,particle(i)%x,particle(i)%y)


        call img2d_norm(a0,FFTsize)
        call img2d_mask(a0,FFTsize,imgmask,.0,.0)
        pr=0.0
        do j=1,FFTsize**2
            pr=pr+a0(j)
        end do
            particle(i)%PR=pr/float(FFTsize**2)
        write(21,900) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
                      particle(i)%x,particle(i)%y,particle(i)%PR,particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1, particle(i)%df2,   particle(i)%astigang!,cent60(3,1),cent60(3,2),cent60(3,3),&
                      !cent60(3,4),cent60(3,5)
    end do
    deallocate(a0)

    close(21)
900 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4)!,f10.4,f10.4,f10.4,f10.4,f10.4
    return
end
