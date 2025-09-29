subroutine ort_cent_13(theta0,phi0,omega0,x0,y0,z0,dang,dtran,pixelsize,ortcent)
    implicit none
    integer jj
    real theta0,phi0,omega0,x0,y0,z0,dang,dtran,ortcent(6,0:12),theta,phi,omega,pixelsize

    ortcent(1,:)=theta0
    ortcent(2,:)=phi0
    ortcent(3,:)=omega0
    ortcent(4,:)=x0
    ortcent(5,:)=y0
    ortcent(6,:)=z0

    ortcent(4,1)=x0-dtran
    ortcent(4,2)=x0+dtran
    ortcent(5,3)=y0-dtran
    ortcent(5,4)=y0+dtran
    ortcent(6,5)=z0-dtran*pixelsize
    ortcent(6,6)=z0+dtran*pixelsize
    ortcent(1,7)=theta0-dang
    ortcent(1,8)=theta0+dang
    ortcent(2,9)=phi0-dang
    ortcent(2,10)=phi0+dang
    ortcent(3,11)=omega0-dang
    ortcent(3,12)=omega0+dang
end

