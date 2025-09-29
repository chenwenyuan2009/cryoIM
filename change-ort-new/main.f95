! A fortran95 program for G95
! By WQY
program main
  implicit none
  integer particleID,imgID,imgparticleID
  integer iargc,n,k1,i,j,sum1
  real    theta,phi,omega,x,y,PR,df1,df2,astigang,dtheta,dphi,domega,x,dy,dPR
  character(len=54):: stckfile
  character(len=128):: ortcentfile
  character(len=64):: par(50)
   n=iargc()
  if (n.lt.2) then
      call help
      goto 999
  end if
  do i=1,n
        call getarg(i,par(i))
  end do
  ortcentfile=trim(par(1))
  mode=int(trim(par(2)))

  open(11,file=trim(ortcentfile))
  if(mode.eq.0)then
  sum1=0
  do i=3,n
     open(21,file=trim(par(i)),status='old')
110  read(21,100,end=120)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     sum1=sum1+1
     write(11,100)stckfile,sum1,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
     goto 110
120  close(21)
    end do
    elseif(mode.eq.1)then
    open(21,file=trim(par(3)),status='old')
130 read(21,100,end=140)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
    write(11,100)stckfile,particleID,theta,phi,omega,x/2.0,y/2.0,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
    goto 130
140 close(21)
    elseif(mode.eq.2)then
    open(21,file=trim(par(3)),status='old')
150 read(21,100,end=160)stckfile,particleID,theta,phi,omega,x,y,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
    write(11,100)stckfile,particleID,theta,phi,omega,x*2.0,y*2.0,PR,imgID,imgparticleID,df1,df2,astigang,&
                         dtheta,dphi,domega,dx,dy,dPR
    goto 150
160 close(21)
    end if
    close(11)
100 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)



subroutine help
    print*,'use change-ort-new newortfile mode=0 ort1name ort2name....... '
    print*,'if mode=1 get half cent :ues change-ort-new newortfile mode=1 oldname'
    print*,'if mode=2 get double cent : use change-ort-new newortfile mode=1 oldname'
end subroutine
999end main
