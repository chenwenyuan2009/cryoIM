subroutine img2d_fill_lin(nlin,a,b,nx,minR,maxR,alpha,FT_lin)
    implicit none
    integer nlin,nx,minR,maxR
    real a(nx**2),b(nx**2),alpha(nlin)
    real FT_lin(2,minR:maxR,nlin)
    real x,y,ax,ay,deg2rad
    integer i,FR,lx,ly,n1,n2

    FT_lin=.0
    deg2rad=atan2(1.0,1.0)/45.
    do i=1,nlin
        if(alpha(i).lt.999.0) then
            do FR=minR,maxR
                x=float(FR)*cos(alpha(i)*deg2rad)+nx/2
                y=float(FR)*sin(alpha(i)*deg2rad)+nx/2

                lx=int(x)
                ax=x-lx
                ly=int(y)
                ay=y-ly
                n1=nx* ly+   lx+1
                n2=nx*(ly+1)+lx+1

                FT_lin(1,FR,i)=(1.0-ax)*(1.0-ay)*a(n1  )+ &
                                    ax *(1.0-ay)*a(n1+1)+ &
                               (1.0-ax)*     ay *a(n2  )+ &
                                    ax *     ay *a(n2+1)
                FT_lin(2,FR,i)=(1.0-ax)*(1.0-ay)*b(n1  )+ &
                                    ax *(1.0-ay)*b(n1+1)+ &
                               (1.0-ax)*     ay *b(n2  )+ &
                                    ax *     ay *b(n2+1)
            end do
        end if
    end do
    return
end subroutine


