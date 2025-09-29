
subroutine FT_amp_expfit(a,b,nx,minR,maxR,Fcut)
    implicit none
    integer nx,minR,maxR,hx,ky,ir,i,Num(minR:maxR)
    real rad,amp
    real a(nx**2),b(nx**2),Fcut(minR:maxR),Favg(minR:maxR),resLN(minR:maxR)
    real aa(2,2),ax(maxR-minR+1,2),bb(2),by(maxR-minR+1)
    real x1,x2

    Favg=0
    Num=0
    resLN=0
    do ky=0,maxR
        do hx=-maxR,maxR
            rad=sqrt(float(hx**2+ky**2))
            ir=int(rad+0.5)
            if((ir.ge.minR).and.(ir.le.maxR)) then
                i=(ky+nx/2)*nx+(hx+nx/2)+1
                amp=sqrt(a(i)**2+b(i)**2)
                Favg(ir)=Favg(ir)+amp
                Num(ir)=Num(ir)+1
                resLN(ir)=resLN(ir)+rad
            end if
        end do
    end do

    do ir=minR,maxR
        if((Num(ir).gt.0).and.(Favg(ir).gt..0)) then
            resLN(ir)=resLN(ir)/float(Num(ir))
            Favg(ir)=alog(Favg(ir)/float(Num(ir)))
        else
            Favg(ir)=.0
        end if
    end do

    ax=1.0
    ax(:,1)=resLN(:)
    by(:)=Favg(:)
    aa=matmul(transpose(ax),ax)
    bb=matmul(transpose(ax),by)
    x2=(bb(2)-bb(1)*aa(2,1)/aa(1,1))/(aa(2,2)-aa(1,2)*aa(2,1)/aa(1,1))
    x1=(bb(1)-aa(1,2)*x2)/aa(1,1)
    do ir=minR,maxR
        Fcut(ir)=exp(resLN(ir)*x1+x2)*0.8
    end do
    return
end


