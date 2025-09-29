subroutine tongji_stck(temp1,nx)
    use common_ort
    integer nx,bin
    character(len=64):: temp1,temp2

    open(31,file=trim(par(4)),status='old')
    bin=0
    nx=0
111    read(31,'(a64)',end=121) temp2
       if(trim(temp2).eq.trim(temp1))then
        bin=bin+1
        goto 111
        endif
        goto 111
121    nx=bin
       close(31)
end




