subroutine img2d_fill_lin_omegearray(Nlin,Nomega,step_angle,a,b,nx,minR,maxR,alpha,FT_lin_omegaarray)
    implicit none
    integer nlin,nx,minR,maxR,Nomega,jj
    real step_angle,a(nx**2),b(nx**2),alpha(nlin)
    real FT_lin_omegaarray(2,minR:maxR,nlin,0:Nomega-1)
    real x,y,ax,ay,deg2rad,omega,tmp,coss,sinn
    integer i,FR,lx,ly,n1,n2

    deg2rad=atan2(1.0,1.0)/45.
    FT_lin_omegaarray=999.0
    do i=1,nlin
        if(alpha(i).lt.999.0) then
            do jj=0,Nomega-1
                omega=alpha(i)-step_angle*jj
                tmp=omega*deg2rad
                coss=cos(tmp)
                sinn=sin(tmp)
                do FR=minR,maxR
                    x=float(FR)*coss+nx/2
                    y=float(FR)*sinn+nx/2

                    lx=int(x)
                    ax=x-lx
                    ly=int(y)
                    ay=y-ly
                    n1=nx* ly+   lx+1
                    n2=nx*(ly+1)+lx+1

                    FT_lin_omegaarray(1,FR,i,jj)=(1.0-ax)*(1.0-ay)*a(n1  )+ &
                                                      ax *(1.0-ay)*a(n1+1)+ &
                                                 (1.0-ax)*     ay *a(n2  )+ &
                                                      ax *     ay *a(n2+1)
                    FT_lin_omegaarray(2,FR,i,jj)=(1.0-ax)*(1.0-ay)*b(n1  )+ &
                                                      ax *(1.0-ay)*b(n1+1)+ &
                                                 (1.0-ax)*     ay *b(n2  )+ &
                                                      ax *     ay *b(n2+1)
                end do
            end do
        end if
    end do
    return
end subroutine


