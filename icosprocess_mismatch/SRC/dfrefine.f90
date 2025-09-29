subroutine dfrefine(df10,df20,astigang0,dtran,dang,dfcent)
    implicit none
    integer jj
    real df10,df20,astigang0,dtran,dang,dfcent(3,7)

    !call return_icos_asymmetry_unit(theta0,phi0,omega0)
    do jj=1,7
        if(jj.eq.1) then
            dfcent(1,jj)=df10
            dfcent(2,jj)=df20
            dfcent(3,jj)=astigang0
        end if
        if(jj.eq.2) then
            dfcent(1,jj)=df10-dtran
            dfcent(2,jj)=df20
            dfcent(3,jj)=astigang0
        end if
        if(jj.eq.3) then
            dfcent(1,jj)=df10+dtran
            dfcent(2,jj)=df20
            dfcent(3,jj)=astigang0
        end if
        if(jj.eq.4) then
            dfcent(1,jj)=df10
            dfcent(2,jj)=df20-dtran
            dfcent(3,jj)=astigang0
        end if
        if(jj.eq.5) then
            dfcent(1,jj)=df10
            dfcent(2,jj)=df20+dtran
            dfcent(3,jj)=astigang0
        end if
        if(jj.eq.6) then
            dfcent(1,jj)=df10
            dfcent(2,jj)=df20
            dfcent(3,jj)=astigang0-dang
        end if
        if(jj.eq.7) then
            dfcent(1,jj)=df10
            dfcent(2,jj)=df20
            dfcent(3,jj)=astigang0+dang
        end if
    end do
    return
end
