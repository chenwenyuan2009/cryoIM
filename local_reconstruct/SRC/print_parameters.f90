subroutine print_parameters
use common_local_reconstruct
implicit none
    print*,'ortfile:       ',trim(ortfile)
    print*,'result3d:      ',trim(result3d)
    print*,'sym:           ', sym
    print*,'Nsym:          ', Nsym
    print*,'FFTsize0:      ', FFTsize0
    print*,'local_FFTsize: ', FFTsize-8
    print*,'totalparticles:', totalparticles
    print*,'useparticle:   ',useparticle
    print*,'apix:          ',apix
    print*,'maxRes:        ',maxRes
    print*,'maxR:          ',maxR
    print*,'imgmask:       ',imgmask
    print*,'last:          ',last
    return
end

