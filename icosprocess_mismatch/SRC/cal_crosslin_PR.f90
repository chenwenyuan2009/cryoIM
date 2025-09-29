! calculate the phase residual of cross-common lines between a particle image and template images according to the cross-common lines of icosahedron.
subroutine cal_crosslin_PR(PR)
    use common_icosproc
    implicit none
    integer i,j,FR
    real theta1,phi1,omega1,theta2,phi2,omega2
    real pi,Tpi,Fcut1(minR:maxR)
    real alpha1(Nlin),alpha2(Nlin)
    real a1(FFTsize**2),b1(FFTsize**2)
    real a2(FFTsize**2),b2(FFTsize**2)
    real FT_commlin1(2,minR:maxR,Nlin),FT_commlin2(2,minR:maxR,Nlin)
    real PR,PR1,PR2,phas1,phas2,amp1,amp2,diff

    pi=4.0*atan2(1.0,1.0)
    Tpi=2.0*pi
    PR=.0
    PR1=.0
    PR2=.0
    a1(:)=template2d(:,1,0)
    b1(:)=template2d(:,2,0)
    call img2d_FT_centshift(a1,b1,FFTsize,template(0)%x,template(0)%y)
    do i=1,totaltemplate
        theta1=template(0)%theta  !image
        phi1  =template(0)%phi
        omega1=template(0)%omega
        theta2=template(i)%theta  !template
        phi2  =template(i)%phi
        omega2=template(i)%omega
        call icos_crosslin_angle(theta1,phi1,omega1,theta2,phi2,omega2,alpha1,alpha2)

        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            a2(:)=template2d(:,1,i)*ctf2d
            b2(:)=template2d(:,2,i)*ctf2d
            call FT_amp_expfit(a2,b2,FFTsize,minR,maxR,Fcut1)
            Fcut(:,i)=Fcut1(:)
        else
            a2(:)=template2d(:,1,i)
            b2(:)=template2d(:,2,i)
        end if
        call img2d_fill_lin(Nlin,a1,b1,FFTsize,minR,maxR,alpha1,FT_commlin1)
        call img2d_fill_lin(Nlin,a2,b2,FFTsize,minR,maxR,alpha2,FT_commlin2)

        do j=1,nlin
            if((alpha1(j).lt.999.).and.(alpha2(j).lt.999.)) then
                do FR=minR,maxR
                    amp1=sqrt(FT_commlin1(1,FR,j)**2+FT_commlin1(2,FR,j)**2)
                    amp2=sqrt(FT_commlin2(1,FR,j)**2+FT_commlin2(2,FR,j)**2)
                    if((amp1.ge.Fcut(FR,0).and.(amp2.ge.Fcut(FR,i)))) then
                        phas1=atan2(FT_commlin1(2,FR,j),FT_commlin1(1,FR,j))
                        phas2=atan2(FT_commlin2(2,FR,j),FT_commlin2(1,FR,j))
                        diff=abs(phas2-phas1)
                        diff=amod(diff,Tpi)
                        if(diff.gt.pi) diff=Tpi-diff
                        amp1=0.5*(amp1+amp2)
                        PR1=PR1+diff*amp1
                        PR2=PR2+amp1
                    end if
                end do
            end if
        end do
    end do

    if(PR2.gt.1.0e-10) then
        PR=PR1/PR2*180.0/pi
    else
        PR=1000.
    end if
    return
end
