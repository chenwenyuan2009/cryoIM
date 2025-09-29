

subroutine img2d_subproj(a1,a2,nx,innerR,outerR)
    implicit none
    integer nx
    real innerR,outerR,fitb0,fitb1,a1(nx**2),a2(nx**2)

    call img2d_fitproj(a1,a2,nx,innerR,outerR,fitb0,fitb1)
    a1=a1*fitb1+fitb0
    a1=a2-a1


    return
 end subroutine







