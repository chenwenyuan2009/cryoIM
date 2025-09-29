subroutine img2d_subproj_new(a1,a2,a3,nx,innerR,outerR)
    implicit none
    integer nx
    real innerR,outerR,fitb0,fitb1,a1(nx**2),a2(nx**2),a3(nx**2)

    call img2d_fitproj(a1,a2,nx,innerR,outerR,fitb0,fitb1)
    a3=a3*fitb1+fitb0
    a3=a2-a3


    return
 end subroutine

