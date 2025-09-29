subroutine get_randort_fromicosort
use common_icosproc
    implicit none
    integer Nparticle,i,n
    real matoular(3,3),maticos(3,3,60),mat(3,3)
    real theta,phi,omega,pi
    real,allocatable::randx(:)

    Nparticle=last-first+1
    allocate(randx(Nparticle))
    call randreal(Nparticle,randx)
    call icos2f_sym_mat(maticos)

    pi=4.0*atan2(1.0,1.0)
    open(11,file='icos_rand.dat')
    do i=first,last
        theta=particle(i)%theta
        phi  =-particle(i)%phi
        omega=180.0-particle(i)%omega
        call euler_mat(theta,phi,omega,matoular)

        n=int(randx(i)*59.0)+1
        if(n.ge.60) n=60
        mat=matmul(maticos(:,:,n),matoular)

        theta=acos(mat(3,3))*180/pi
        phi  =-atan2(mat(2,3),mat(1,3))*180/pi
        omega=atan2(mat(3,2),mat(3,1))*180/pi

        particle(i)%theta=theta
        particle(i)%phi  =phi
        particle(i)%omega=omega

        write(11,100) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
        particle(i)%x,particle(i)%y,particle(i)%PR, particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1,       particle(i)%df2,   particle(i)%astigang, &
                      particle(i)%dtheta,particle(i)%dphi,particle(i)%domega, &
                      particle(i)%dx,particle(i)%dy,particle(i)%dPR
    end do
100 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
    close(11)
end

subroutine randreal(n,randx)
    implicit none
    integer n,i
    real randx(n)
    randx(:)=.0
    call init_random_seed_use_timer()
    do i=1,n
        call random_number(randx(i))
    end do
end

subroutine init_random_seed_use_timer()
        integer :: i, n, clock
        integer, allocatable :: seed(:)

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)

        deallocate(seed)
    end subroutine






