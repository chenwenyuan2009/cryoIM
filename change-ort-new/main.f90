program main
  use common_ort
  implicit none
  integer particleID,imgID,imgparticleID,IDp,local_x,local_y
  integer iargc,n,k1,i,j,sum1,Nmax,k,jj,first,last,FFTsize,bin1
  real    theta,phi,omega,x,y,PR,df1,df2,astigang,dtheta,dphi,domega,dx,dy,dPR
  real    theta1,phi1,omega1,centx,centy,local_z,PR1,alpha,phi2,bound,pr2
  real    x0,y0,z0,initiald(3,1),RotMat(3,3),newd(3,1),ix0,iy0,radius
  character(len=64) stckfile2,suffix


   n=iargc()
  if (n.lt.3) then
      call help
      goto 999
  end if
  do i=1,n
        call getarg(i,par(i))
  end do
  ortcentfile=trim(par(1))
  call input_parameters
  do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'mode') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) mode
                end if
                 if(trim(coeff(j)).eq.'FFTsize') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) FFTsize
                end if
            end if

        end do
   end do
  print*,'mode=',mode

  open(11,file=trim(ortcentfile))
   sum1=0
  if(mode.eq.0)then

     do i=3,n
     open(21,file=trim(par(i)),status='old')
110  read(21,*,end=120)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     sum1=sum1+1

     write(11,100)stckfile,sum1,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 110
120  close(21)
     print*,sum1
     end do
    elseif(mode.eq.6)then
    sum1=0
     do i=3,n
     open(21,file=trim(par(i)),status='old')
111  read(21,*,end=121)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang

     sum1=sum1+1

     write(11,101)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang

     goto 111
121  close(21)
     print*,sum1
     end do
   elseif(mode.eq.1)then
     open(21,file=trim(par(3)),status='old')
130  read(21,*,end=140)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
                         sum1=sum1+1
                         print*,sum1
     write(11,100)stckfile,particleID,theta,phi,omega,(x-0.5)/2.0,(y-0.5)/2.0,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 130
140  close(21)
   elseif(mode.eq.2)then
     open(21,file=trim(par(3)),status='old')
150  read(21,*,end=160)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
                        sum1=sum1+1
                         print*,sum1
        k=index(stckfile,char(46),BACK=.TRUE.)
        temp=stckfile(1:k-1)//'.mrcs'
     write(11,100)temp,particleID,theta,phi,omega,x*2.0+0.5,y*2.0+0.5,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 150
160  close(21)
   elseif(mode.eq.8)then
     open(21,file=trim(par(3)),status='old')
250  read(21,*,end=260)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
        k=index(stckfile,char(46),BACK=.TRUE.)
        temp=stckfile(1:k-1)//'.bin2'
     write(11,100)temp,particleID,theta,phi,omega,x*2.0+0.5,y*2.0+0.5,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 250
260  close(21)
    elseif(mode.eq.3)then
     open(21,file=trim(par(3)),status='old')
170  read(21,*,end=180)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     write(11,100)stckfile,particleID,theta,phi,omega,x,y,cos(PR*3.1415926/180.0),imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 170
180  close(21)
     elseif(mode.eq.4)then
     open(21,file=trim(par(3)),status='old')
190  read(21,*,end=200)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,dtheta
     write(11,101)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,(df1+df2)/2.0,(df1+df2)/2.0,0.0,0.1

     goto 190
200  close(21)
     elseif(mode.eq.5)then
     open(21,file=trim(par(3)),status='old')
     write(*,'(a6)'), 'bound','PR'
     read(*,*) bound,pr2
210  read(21,*,end=220)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang!,&
                         !dtheta,dphi,domega,dx,dy,dPR
                         if(abs(x).lt.bound .and. abs(y).lt. bound.and. pr.lt.pr2)then
                            sum1=sum1+1
     write(11,100)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
!                         dtheta,dphi,domega,dx,dy,dPR
                         end if

     goto 210
220  close(21)
     print*,sum1
     elseif(mode.eq.7)then
     open(21,file=trim(par(3)),status='old')
230  read(21,*,end=240)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang!,&
                         !dtheta,dphi,domega,dx,dy,dPR
     write(11,102)      particleID,theta,-phi,omega,x,y,PR,imgID,df1,df2,astigang
     goto 230
240  close(21)
     elseif(mode.eq.9)then
     open(21,file=trim(par(3)),status='old')
270  read(21,*,end=280)particleID,theta,phi,omega,x,y,PR,imgID,df1,df2,astigang
     write(11,102)  particleID,theta,-phi,omega,x,y,PR,imgID,df1,df2,astigang
     goto 270
280  close(21)
     elseif(mode.eq.10)then
     open(21,file=trim(par(3)),status='old')
290  read(21,*,end=300)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
     write(11,101)  stckfile,particleID,theta,phi,omega,-y,-x,PR,imgID,imgparticleID,df1,df2,astigang
     goto 290
300  close(21)
     elseif(mode.eq.11)then
     open(21,file=trim(par(3)),status='old')
310  read(21,*,end=320)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
     write(11,101)  stckfile,particleID,theta,phi,omega,x,y,1.0-PR,imgID,imgparticleID,df1,df2,astigang
     goto 310
320  close(21)
     elseif(mode.eq.12)then
        stckfile1='abcdef'
     open(21,file=trim(par(3)),status='old')
330  read(21,*,end=340)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
     if(stckfile1.ne.trim(stckfile2))then
     stckfile1=trim(stckfile2)
     call tongji_stck(stckfile1,jj)

     bin1=0
     endif
     bin1=bin1+1
     write(11,1010)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,bin1,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1,jj
     goto 330
340  close(21)
     elseif(mode.eq.13)then
     open(21,file=trim(par(3)),status='old')
     write(*,*) 'please input x0 y0 z0 (pixel) radius'
     read(*,*) x0,y0,z0,radius
     initiald(1,1)=x0;initiald(2,1)=y0;initiald(3,1)=z0
350  read(21,*,end=360)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang
     if(abs(x).le. radius .and. abs(y) .le. radius)then
     call oular_mat(RotMat,theta,phi,omega)
     newd=matmul(RotMat,initiald)
     ix0=newd(1,1)+x
     iy0=newd(2,1)+y


     write(11,103) stckfile,particleID,theta,phi,omega,ix0-int(ix0),iy0-int(iy0),PR,imgID,imgparticleID&
                   ,df1,df2,astigang,ix0,iy0
     endif
     goto 350
360  close(21)

     elseif(mode.eq.14)then
     open(21,file=trim(par(3)),status='old')
370  read(21,*,end=380)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
     IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
     write(11,104)  stckfile2,particleID,theta,phi,omega,(x-0.5)/2.0,(y-0.5)/2.0,PR,imgID,imgparticleID,&
     df1,df2,astigang,&
     IDp,(local_x+1)/2,(local_y+1)/2,theta1,phi1,omega1,(centx-0.5)/2.0,(centy-0.5)/2.0,local_z,PR1
     goto 370
380  close(21)
     elseif(mode.eq.15)then
     open(21,file=trim(par(3)),status='old')
     print*,'please input suffix_name such as .bin2'
     read(*,*) suffix
390  read(21,*,end=400)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
     IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
        k=index(stckfile2,char(46),BACK=.TRUE.)
        temp=stckfile2(1:k-1)//trim(suffix)
     write(11,104)  temp,particleID,theta,phi,omega,x*2+0.5,y*2+0.5,PR,imgID,imgparticleID,df1,df2,&
     astigang,IDp,(local_x)*2-1,(local_y)*2-1,theta1,phi1,omega1,centx*2+1.0,centy*2.0+1.0,local_z,PR1
     goto 390
400  close(21)

     elseif(mode.eq.16)then
     open(21,file=trim(par(3)),status='old')
     !open(22,file=trim(par(4)),status='old')
     write(*,*)'please input copy_number'
     read(*,*) first
     alpha=360.0/float(first)

 303      read(21,*,end=302)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       do j=1,first
        phi2=(j-1)*alpha+phi1
        if(phi2.gt.360.0) phi2=phi2-360.0
       write(11,104)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi2,omega1,centx,centy,local_z,PR1
       end do
       goto 303
302     close(21)
     !close(22)

     elseif(mode.eq.17)then
        print*,'please input PR_threshold'
        read(*,*) x0
     open(21,file=trim(par(3)),status='old')
     open(22,file=trim(par(4)))
301    read(21,*,end=311)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       if(PR.lt.x0)then
       write(11,104)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 301
       else
       write(22,104)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 301
       endif

311  close(21)
     close(22)

    elseif(mode.eq.18)then
        print*,'please input PR_threshold'
        read(*,*) x0
     open(21,file=trim(par(3)),status='old')
     open(22,file=trim(par(4)))
410    read(21,*,end=411)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang

       if(PR.lt.x0)then
       write(11,101)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang

       goto 410
       else
       write(22,101)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang

       goto 410
       endif

411  close(21)
     close(22)

     elseif(mode.eq.19)then

     open(21,file=trim(par(3)),status='old')
432    read(21,*,end=431)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       local_x=local_x+int(centx)
       local_y=local_y+int(centy)
       centx=centx-int(centx)
       centy=centy-int(centy)
       write(11,104)  stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 432

431  close(21)
     !close(22)

     elseif(mode.eq.20)then
     open(21,file=trim(par(3)),status='old')
442    read(21,*,end=441)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       x=float(local_x)+centx-FFTsize/2.0
       y=float(local_y)+centy-FFTsize/2.0
       df1=df1/10000.0
       df2=df2/10000.0
       write(11,104)  stckfile2,particleID,theta1,phi1,omega1,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 442

441  close(21)
     !close(22)
     elseif(mode.eq.21)then
     open(21,file=trim(par(3)),status='old')
     open(22,file=trim(par(4)),status='old')
     stckfile='abcdef'
452    read(21,*,end=451)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1

       if(trim(stckfile2).ne.trim(stckfile))then
199      read(22,*,end=81) stckfile1
          if(trim(stckfile1).eq.trim(stckfile2))then
            stckfile=trim(stckfile1)
       write(11,104)  stckfile2,particleID,theta1,phi1,omega1,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       rewind(22)
 81    goto 452
      else
        goto 199
        endif
      else
      write(11,104)  stckfile2,particleID,theta1,phi1,omega1,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
     IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 452
451  close(21)

     endif
     close(21)
     close(22)

     elseif(mode.eq.22)then
        print*, 'please input 0 or 1, 0 for 除以一万，1 for 乘以一万'
        read(*,*) i
     open(21,file=trim(par(3)),status='old')
453    read(21,*,end=454)stckfile2,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       if(i.eq. 0)then
       df1=df1/10000.0
       df2=df2/10000.0
       elseif(i.eq. 1)then
        df1=df1*10000.0
        df2=df2*10000.0
        endif
       write(11,104)  stckfile2,particleID,theta1,phi1,omega1,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
       IDp,local_x,local_y,theta1,phi1,omega1,centx,centy,local_z,PR1
       goto 453
454  close(21)
     endif
     close(11)
102 format(i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,f15.4,f15.4,f10.4)
101 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f15.4,f15.4,f10.4)
103 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f15.4,f15.4,f10.4,f10.4,f10.4)
1010 format(a64,I9,F10.4,F10.4,F10.4,  &
           F10.4,F10.4,F10.4,I11,I7,  &
           F15.4,F15.4,F10.4, &
           i9,i6,i6,  &
           f10.4,f10.4,f10.4,  &
           f10.4,f10.4,f11.4,f10.4,I4)
100 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f15.4,f15.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
104 format(a64,I9,F10.4,F10.4,F10.4,  &
           F10.4,F10.4,F10.4,I11,I7,  &
           F15.4,F15.4,F10.4, &
           i9,i6,i6,  &
           f10.4,f10.4,f10.4,  &
           f10.4,f10.4,f11.4,f10.4)

999 end

