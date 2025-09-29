subroutine input_parameters
use common_local_reconstruct
    coeff(1) ='local_ort'
    coeff(2) ='local_sym'
    coeff(3) ='local_FFTsize'
    coeff(4) ='apix'
    coeff(5) ='maxRes'
    coeff(6) ='imgmask'  !angstrom
    coeff(7) ='first'
    coeff(8) ='last'
    coeff(9) ='result3d'
    coeff(10) ='applyCTF'
    coeff(11) ='vol'
    coeff(12) ='Cs'
    coeff(13)='Bfactor'
    coeff(14)='ampwg'
    coeff(15)='PR_threshold'
    coeff(16)='subrec'
    coeff(17)='imgstck'
    coeff(18)='boundX'
    coeff(19)='Zhigh_df'
    coeff(20)='sign'
    coeff(21)='tilt'

    proc2d%applyCTF='y'
    proc2d%imgmask ='n'
    proc2d%norm    ='y'
    proc2d%imgstck ='n'
    proc2d%last    ='n'
    proc3d%subrec  ='n'
    proc3d%Zhigh_df='n'
    local_sym='c1'
    result3d='local_map.mrc'
    PR_threshold  =90.0
    CTFpar%vol    =200
    CTFpar%Cs     =2.7
    CTFpar%ampwg  =0.1
    CTFpar%Bfactor=.0
    first         =1
    Nsym=1
    FFTsize=256  !local FFTsize
    boundX =100
    maxRes=10.0
    apix=1.0
    tilt=90.0
    return
end
