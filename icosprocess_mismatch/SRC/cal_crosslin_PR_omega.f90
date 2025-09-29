

subroutine cal_crosslin_PR_omega(Nomega,PR1_omega)
    use common_icosproc
    implicit none
    integer i,j,FR,Nomega,iomega
    real theta1,phi1,omega1,theta2,phi2,omega2
    real rad2deg,pi,Tpi,phas1,phas2,amp1,amp2,diff
    real alpha1(Nlin),alpha2(Nlin),Fcut1(minR:maxR)
    real a1(FFTsize**2),b1(FFTsize**2)
    real a2(FFTsize**2),b2(FFTsize**2)
    real,allocatable::FT_commlin1_omega(:,:,:,:),FT_commlin2(:,:,:)
    real PR1_omega(0:Nomega-1),PR2_omega(0:Nomega-1)
!    allocate(FT_commlin1(2,minR:maxR,Nlin))
    allocate(FT_commlin1_omega(2,minR:maxR,nlin,0:Nomega))
    allocate(FT_commlin2(2,minR:maxR,Nlin))

    pi =4.0*atan2(1.0,1.0)
    Tpi=8.0*atan2(1.0,1.0)
    rad2deg=180.0/pi
    PR1_omega=.0
    PR2_omega=.0
    FT_commlin1_omega=.0
    a1(:)=template2d(:,1,0)
    b1(:)=template2d(:,2,0)
!    call img2d_FT_centshift(a1,b1,FFTsize,template(0)%x,template(0)%y)
    do i=1,totaltemplate
        theta1=template(0)%theta  !image
        phi1  =template(0)%phi
        omega1=template(0)%omega
        theta2=template(i)%theta  !template
        phi2  =template(i)%phi
        omega2=template(i)%omega
        call icos_crosslin_angle(theta1,phi1,omega1,theta2,phi2,omega2,alpha1,alpha2) !calculate the angle according to cross-correlation
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            a2(:)=template2d(:,1,i)*ctf2d
            b2(:)=template2d(:,2,i)*ctf2d
            call FT_amp_expfit(a2,b2,FFTsize,minR,maxR,Fcut1)
            Fcut(:,i)=Fcut1(:)
        else
            a2(:)=template2d(:,1,i)
            b2(:)=template2d(:,2,i)
        end if
        call img2d_fill_lin(Nlin,a2,b2,FFTsize,minR,maxR,alpha2,FT_commlin2)
        call img2d_fill_lin_omegearray(Nlin,Nomega,searchstep,a1,b1,FFTsize,minR,maxR,alpha1,FT_commlin1_omega)

        !!!!
!            call img2d_fill_lin(Nlin,a1,b1,FFTsize,minR,maxR,alpha1-step_angle*iomega,FT_commlin1)
        do j=1,nlin
            if((alpha1(j).lt.999.).and.(alpha2(j).lt.999.)) then
                do FR=minR,maxR
                    amp2=sqrt(FT_commlin2(1,FR,j)**2+FT_commlin2(2,FR,j)**2)
                    if(amp2.ge.Fcut(FR,i)) then
                        phas2=atan2(FT_commlin2(2,FR,j),FT_commlin2(1,FR,j))
                        do iomega=0,Nomega-1
                            amp1=sqrt(FT_commlin1_omega(1,FR,j,iomega)**2+FT_commlin1_omega(2,FR,j,iomega)**2)
                            if(amp1.ge.Fcut(FR,0)) then
                                phas1=atan2(FT_commlin1_omega(2,FR,j,iomega),FT_commlin1_omega(1,FR,j,iomega))
                                diff=abs(phas2-phas1)
                                diff=amod(diff,Tpi)
                                if(diff.gt.pi) diff=Tpi-diff
                                amp1=0.5*(amp1+amp2)
                                PR1_omega(iomega)=PR1_omega(iomega)+diff*amp1
                                PR2_omega(iomega)=PR2_omega(iomega)+amp1
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end do

    do iomega=0,Nomega-1
        if(PR2_omega(iomega).gt.1.0e-10) then
            PR1_omega(iomega)=PR1_omega(iomega)/PR2_omega(iomega)*rad2deg
        else
            PR1_omega(iomega)=999.0
        end if
    end do
    deallocate(FT_commlin1_omega)
    deallocate(FT_commlin2)
    return
end
