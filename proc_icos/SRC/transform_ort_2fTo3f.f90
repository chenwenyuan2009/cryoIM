subroutine transform_ort_2fTo3f
use common_proc_icos
implicit none
    integer particleID,imgID,imgparticleID,i1
    real theta,phi,omega,r(3,3),v,u,w,deg2rad,rot0(3,3),rad2deg

    open(12,file=trim(newort))

    deg2rad=4.0*atan2(1.0,1.0)/180.0
    rad2deg=180.0/4.0/atan2(1.0,1.0)
    v=acos(cos(60.0*deg2rad)/sin(36.0*deg2rad))     !angle between 2f and 5f
    u=acos(1.0/tan(60.0*deg2rad)/tan(36.0*deg2rad)) !angle between 3f and 5f
    w=acos(cos(36.0*deg2rad)/sin(60.0*deg2rad))     !angle between 2f and 3f
    rot0(1,1)=1.; rot0(1,2)=.0     ; rot0(1,3)=.0
    rot0(2,1)=.0; rot0(2,2)= cos(-w); rot0(2,3)=sin(-w)
    rot0(3,1)=.0; rot0(3,2)=-sin(-w); rot0(3,3)=cos(-w)
    call load_original_ortcent 
    do i1=1,totalparticle
        call Eularmatrix(particle(i1)%theta,particle(i1)%phi,particle(i1)%omega,r)
        r=matmul(rot0,r)
        if(r(3,3).ge.1.0) then
            theta=0
        else if(r(3,3).lt.-1.0) then
            theta=180.0
        else
            theta=acos(r(3,3))*rad2deg
        end if
        phi  =atan2(r(2,3), r(1,3))*rad2deg
        omega=atan2(r(3,2),-r(3,1))*rad2deg
        particle(i1)%theta=theta
        particle(i1)%phi  =phi
        particle(i1)%omega=omega
        write(12,900)  particle(i1)%stckfile,particle(i1)%particleID,particle(i1)%theta,particle(i1)%phi,particle(i1)%omega, &
                       particle(i1)%x,particle(i1)%y,particle(i1)%PR, particle(i1)%imgID,particle(i1)%imgparticleID,  &
                       particle(i1)%df1,particle(i1)%df2,particle(i1)%astigang
    end do

900 format(a64,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f15.4,f15.4,f10.4)
998 print*,'the new ortfile is written in ',trim(newort)
999 end 
