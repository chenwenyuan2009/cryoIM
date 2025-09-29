!存储各输入参数

subroutine input_parameters
use common_icos_loc_box
    implicit none
    coeff(1) ='icos_ort'
    coeff(2) ='local_ort'
    coeff(3) ='local_imgstck'
    coeff(4) ='local_x'
    coeff(5) ='local_y'
    coeff(6) ='local_z'
    coeff(7) ='apix'
    coeff(8) ='PR_threshold'
    coeff(9) ='sub_FFTsize'
    coeff(10)='icos_imgstck'
    coeff(11)='Nsym'
    coeff(12)='sym_matrix'
    coeff(13)='sym'
    coeff(14)='first'
    coeff(15)='last'
    coeff(16)='boundX'

    first=1
    apix=1
    PR_threshold=90
    sub_FFTsize=256
    proc2d_icos_imgstck='n'
    proc2d_last='n'
    sym='i2'
    Nsym=60
    sym_matrix='sym_matrix.txt'
    proc2d_sym_matrix='n'
    proc2d_local_imgstck='n'
    boundX=100
    return
end

