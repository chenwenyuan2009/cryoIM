subroutine input_parameters
use common_loc_refine
    implicit none

    coeff(1) ='local_ort0'
    coeff(2) ='local_newort'
    coeff(3) ='model3d'
    coeff(4) ='imgstck'
    coeff(5) ='minRes'
    coeff(6) ='apix'
    coeff(7) ='first'
    coeff(8) ='last'
    coeff(9) ='imgmask' !angstrom
    coeff(10) ='applyCTF'
    coeff(11)='Ncycle'
    coeff(12)='vol'
    coeff(13)='Cs'
    coeff(14)='refine_angle'
    coeff(15)='refine_shift'
    coeff(16)='maxRes'
    coeff(17)='convF'
    coeff(18)='particle_df'
    coeff(19)='Zhigh_df'
    coeff(20)='sym'
    coeff(21)='sign'
    coeff(22)='ampwg'

    proc2d%norm='y'
    proc2d%imgmask='n'
    proc2d%applyCTF='y'
    proc2d%imgstck='n'
    proc2d%centshift='y'
    proc2d%first='n'
    proc2d%last='n'
    proc2d%applyCTF='y'
    proc2d%particle_df='n'
    proc2d%Zhigh_df='n'

    CTFpar%vol   =200
    CTFpar%Cs    =2.7
    CTFpar%ampwg =.1
    CTFpar%Bfactor=.0
    proc2d%refine_angle='n'
    proc2d%refine_shift='n'

    first=1
    apix=1.0
    refine_angle=1.0
    refine_shift=1.0
    Ncycle=12
    minRes=200.0   !angstrom
    maxRes=10.0
    convF=0.8
    sym='c1'
    return
end

