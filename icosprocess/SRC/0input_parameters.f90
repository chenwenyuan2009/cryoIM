!存储各输入参数

subroutine input_parameters
use common_icosproc
    implicit none
    coeff(1) ='apix'
    coeff(2) ='minRes'
    coeff(3) ='maxRes'
    coeff(4) ='imgmask'
    coeff(5) ='first'
    coeff(6) ='last'
    coeff(7) ='applyCTF'
    coeff(8) ='Ncycle'
    coeff(9)='templatestck'
    coeff(10)='vol'
    coeff(11)='Cs'
    coeff(12)='Bfactor'
    coeff(13)='ampwg'
    coeff(14)='PR_threshold'
    coeff(15)='result3d'
    coeff(16)='mode'
    coeff(17)='newortfile'
    coeff(18)='model3d'
    coeff(19)='templateort'
    coeff(20)='searchstep'
    coeff(21)='realspaceavg'
    coeff(22)='ortfile'
    coeff(23)='recISAF'
    coeff(24)='subrec'
    coeff(25)='FSC'
    coeff(26)='imgstck'
    coeff(27)='sym'
    coeff(28)='maxcentshift'
    coeff(29)='corrCTF'
    coeff(30)='boundX'
    coeff(31)='loc_box'
    coeff(32)='check'


    proc2d%applyCTF    ='y'
    proc2d%check       ='n'
    proc2d%imgmask     ='y'
    proc2d%corrCTF     ='n'
    proc3d%realspaceavg='n'
    proc3d%model3d     ='n'
    proc3d%templatestck='n'
    proc3d%templateort ='n'
    proc3d%recISAF     ='n'
    proc3d%subrec      ='n'
    sym='icos'
    maxcentshift=100
    return
end

