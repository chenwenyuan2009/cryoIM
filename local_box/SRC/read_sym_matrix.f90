subroutine read_sym_matrixa(sym_matrix_file,sym_rot,N)
implicit none
    integer N,i,j
    real sym_rot(3,3,N)
    character(len=64)::sym_matrix_file

    open(11,file=trim(sym_matrix_file),status='old')
    do i=1,N
        read(11,*) sym_rot(1,1,i),sym_rot(2,1,i),sym_rot(3,1,i)
        read(11,*) sym_rot(1,2,i),sym_rot(2,2,i),sym_rot(3,2,i)
        read(11,*) sym_rot(1,3,i),sym_rot(2,3,i),sym_rot(3,3,i)
    end do
    close(11)
    return
    end
