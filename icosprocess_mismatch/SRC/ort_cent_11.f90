subroutine ort_cent_11(theta0,phi0,omega0,x0,y0,dang,dtran,ortcent)
    implicit none
    integer jj
    real theta0,phi0,omega0,x0,y0,dang,dtran,ortcent(5,11),theta,phi,omega

    !call return_icos_asymmetry_unit(theta0,phi0,omega0)
    do jj=1,11
        if(jj.eq.1) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.2) then
            ortcent(1,jj)=theta0-dang
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.3) then
            ortcent(1,jj)=theta0+dang
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.4) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0-dang
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.5) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0+dang
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.6) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0-dang
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.7) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0+dang
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0
        end if
        if(jj.eq.8) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0-dtran
            ortcent(5,jj)=y0
        end if
        if(jj.eq.9) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0+dtran
            ortcent(5,jj)=y0
        end if
        if(jj.eq.10) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0-dtran
        end if
        if(jj.eq.11) then
            ortcent(1,jj)=theta0
            ortcent(2,jj)=phi0
            ortcent(3,jj)=omega0
            ortcent(4,jj)=x0
            ortcent(5,jj)=y0+dtran
        end if
        theta=ortcent(1,jj)
        phi  =ortcent(2,jj)
        omega=ortcent(3,jj)
        !call return_icos_asymmetry_unit(theta,phi,omega)
        ortcent(1,jj)=theta
        ortcent(2,jj)=phi
        ortcent(3,jj)=omega
    end do
    return
end

