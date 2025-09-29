subroutine randomreal(randx,N)
    implicit none
    integer i,N
    real randx(N)
    randx(:)=.0
    call init_random_seed_use_timer()
    do i=1,N
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




