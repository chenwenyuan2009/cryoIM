subroutine load_template
use common_icosproc
    implicit none
    integer i,hx,ky
    character(len=54)::char0

    real,allocatable::shift2d(:),a(:),b(:),Fcut1(:)
    allocate(shift2d(FFTsize**2),a(FFTsize**2),b(FFTsize**2),Fcut1(minR:maxR))

   open(21,file=trim(templateort), status='old')
   print*,'loading 2D templates...'
    do i=1,FFTsize*FFTsize
        ky=(i-1)/FFTsize
        hx=mod(i-1,FFTsize)
        shift2d(i)=(-1)**(hx+ky)
    end do
    
    do i=1,totaltemplate
        read(21,*) template(i)%particleID,template(i)%theta,template(i)%phi,template(i)%omega, &
                   template(i)%x,template(i)%y
    end do
    close(21)

    if(proc3d%templatestck.eq.'y') then
        open(31,file=trim(templatestck),form='unformatted',access='stream',status='old')
        read(31) mrc
        do i=1,totaltemplate
            read(31) template2d(:,1,i)
        end do
        close(31)
    else
        call project_template2d
    end if

    do i=1,totaltemplate
        a=template2d(:,1,i)
        if((proc2d%imgmask.eq.'y').or.(proc2d%imgmask.eq.'Y')) call img2d_mask(a,FFTsize,imgmask,template(i)%x,template(i)%y)
        a=a*shift2d
        b=.0
        call img2d_FFT(a,b,FFTsize,FFTsize,1)
        call img2d_FT_centshift(a,b,FFTsize,template(i)%x,template(i)%y)
        template2d(:,1,i)=a(:)
        template2d(:,2,i)=b(:)
        if((proc2d%applyCTF.ne.'y').or.(proc2d%applyCTF.ne.'Y')) then
            call FT_amp_expfit(a,b,FFTsize,minR,maxR,Fcut1)
            Fcut(:,i)=Fcut1(:)
        end if
    end do
    deallocate(shift2d,a,b,Fcut1)
    return
end subroutine




