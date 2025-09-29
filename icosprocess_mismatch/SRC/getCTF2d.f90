
subroutine getCTF2d(nx,df10,df20,Astigangle0,apix,vol,Cs0,ampwg,CTF2d)
    implicit none
    integer nx,i,hx,ky
    real df10,df20,Astigangle0,Cs0
    real df1,df2,df,Astigangle,apix,vol,Cs,ampwg,lamada
    real deg2rad,PI,phi,res,phase
    real CTF2d(nx**2)

    PI=4.0*atan2(1.0,1.0)
    deg2rad=PI/180.0
    lamada= 0.387818/sqrt(1.0*vol*(1.0+0.978459E-3*vol))
    Cs=Cs0*1.0e7
    df1=df10*1.0e4
    df2=df20*1.0e4
    Astigangle=Astigangle0*deg2rad
!    print*,Astigangle,lamada
!    print*,PI
!    print*,deg2rad
!    print*,Cs0

    do i=1,nx**2
        ky=(i-1)/nx-nx/2
        hx=mod(i-1,nx)-nx/2
        phi=atan2(float(ky),float(hx))
        if(phi.lt..0) phi=phi+2.0*PI
        phi=phi-Astigangle
        df=df1*cos(phi)**2+df2*sin(phi)**2
        res=sqrt(float(hx**2+ky**2))/float(nx)/apix
        phase=2.0*PI*(-cs*lamada**3*res**4/4+df*lamada*res**2/2)
        ctf2d(i)=-(sqrt(1.0-ampwg**2)*sin(phase)+ampwg*cos(phase))
      end do
    return
end


