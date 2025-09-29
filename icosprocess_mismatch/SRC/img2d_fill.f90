subroutine img2d_fill(a,b,nx,minR,maxR,alpha,FT_fill)
    implicit none
    integer nx,minR,maxR
    real a(nx**2),b(nx**2)
    real FT_fill(2,minR:maxR),alpha
    real x,y,ax,ay,deg2rad
    integer FR,lx,ly,n1,n2

    FT_fill=.0
    deg2rad=atan2(1.0,1.0)/45.0
    if(alpha.lt.999.0)then
            do FR=minR,maxR
                x=float(FR)*cos(alpha*deg2rad)+nx/2
                y=float(FR)*sin(alpha*deg2rad)+nx/2

                lx=int(x)
                ax=x-lx
                ly=int(y)
                ay=y-ly
                n1=nx* ly+   lx+1
                n2=nx*(ly+1)+lx+1

                FT_fill(1,FR)=(1.0-ax)*(1.0-ay)*a(n1  )+ &
                                    ax *(1.0-ay)*a(n1+1)+ &
                               (1.0-ax)*     ay *a(n2  )+ &
                                    ax *     ay *a(n2+1)
                FT_fill(2,FR)=(1.0-ax)*(1.0-ay)*b(n1  )+ &
                                    ax *(1.0-ay)*b(n1+1)+ &
                               (1.0-ax)*     ay *b(n2  )+ &
                                    ax *     ay *b(n2+1)
                end do
            end if
    return
end subroutine
