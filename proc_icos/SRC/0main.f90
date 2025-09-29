
program main
    use common_proc_icos
    implicit none
    integer iargc,n,k1,i,j

    n=iargc()
    if (n.lt.1) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
    end do

    call inputpars
    do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'ort0') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ort0
                end if
                if(trim(coeff(j)).eq.'newort') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) newort
                end if
                if(trim(coeff(j)).eq.'calc_randort_2f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%calc_randort_2f
                end if
                if(trim(coeff(j)).eq.'calc_matrix_2f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%calc_matrix_2f
                end if
                if(trim(coeff(j)).eq.'calc_matrix_3f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%calc_matrix_3f
                end if
                if(trim(coeff(j)).eq.'calc_matrix_5f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%calc_matrix_5f
                end if
                if(trim(coeff(j)).eq.'transform_ort_2fTo5f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%transform_ort_2fTo5f
                end if
                if(trim(coeff(j)).eq.'transform_ort_2fTo3f') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc%transform_ort_2fTo3f
                end if

                goto 90
            end if
        end do
        call help
        goto 999
90 end do

    if((proc%calc_matrix_2f.eq.'Y').or.(proc%calc_matrix_2f.eq.'y')) call calc_matrix_2f
    if((proc%calc_matrix_3f.eq.'Y').or.(proc%calc_matrix_3f.eq.'y')) call calc_matrix_3f
    if((proc%calc_matrix_5f.eq.'Y').or.(proc%calc_matrix_5f.eq.'y')) call calc_matrix_5f
    if((proc%calc_randort_2f.eq.'Y').or.(proc%calc_randort_2f.eq.'y')) call calc_randort_2f
    if((proc%transform_ort_2fTo5f.eq.'Y').or.(proc%transform_ort_2fTo5f.eq.'y')) call transform_ort_2fTo5f
    if((proc%transform_ort_2fTo3f.eq.'Y').or.(proc%transform_ort_2fTo3f.eq.'y')) call transform_ort_2fTo3f

999 end


