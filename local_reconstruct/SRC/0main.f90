!This package is used to cryo-EM data process for icosahedral virus.
!The program includes six modules:
!             0.2D project. need input 3D model and orientation files.
!             1.3D reconstruction.
!             2.Random initial model generation.
!             3.Initial particle center determination.
!             4.Icosahedral orientation globe search.
!             5.Icosahedral orientation local refinement.
!Hongrong Liu. 03/15/2017.  hrliu@hunnu.edu.cn
!Intel Fortran compile 03/15/2017

program main
    use common_local_reconstruct
    implicit none
    integer iargc,n,k1,i,j,particleID
    real theta,phi,omega,x0,y0
    character(len=64)::char0


    n=iargc()
    if (n.lt.1) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
    end do

    call input_parameters
    do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'local_ort') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ortfile
                end if
                if(trim(coeff(j)).eq.'local_sym') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym
                end if
                if(trim(coeff(j)).eq.'local_FFTsize') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) FFTsize
                end if
                if(trim(coeff(j)).eq.'apix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) apix
                end if
                if(trim(coeff(j)).eq.'maxRes') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) maxRes
                end if
                if(trim(coeff(j)).eq.'imgmask') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgmask
                    if(imgmask.gt.0) proc2d%imgmask='y'
                end if
                if(trim(coeff(j)).eq.'first') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) first
                end if
                if(trim(coeff(j)).eq.'last') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) last
                    proc2d%last='y'
                end if
                if(trim(coeff(j)).eq.'result3d') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) result3d
                end if
                if(trim(coeff(j)).eq.'applyCTF') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%applyCTF
                end if
                if(trim(coeff(j)).eq.'vol') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) CTFpar%vol
                end if
                if(trim(coeff(j)).eq.'Cs') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) CTFpar%Cs
                end if
                if(trim(coeff(j)).eq.'ampwg') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) CTFpar%ampwg
                end if
                if(trim(coeff(j)).eq.'Bfactor') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) CTFpar%Bfactor
                end if
                if(trim(coeff(j)).eq.'PR_threshold') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) PR_threshold
                end if
                if(trim(coeff(j)).eq.'boundX') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) boundX
                end if
                if(trim(coeff(j)).eq.'subrec') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%subrec
                end if
                if(trim(coeff(j)).eq.'Zhigh_df') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%Zhigh_df
                end if
                if(trim(coeff(j)).eq.'imgstck') then
                    proc2d%imgstck='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgstck
                end if
                if(trim(coeff(j)).eq.'sign') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sign
                end if
                if(trim(coeff(j)).eq.'tilt')then
                    read(par(i)(k1+1:len(trim(par(i)))),*) tilt
                end if
                goto 110
            end if
        end do
        call help
        goto 999
110 end do
    
    call calc_initialvalues
    print*,imgmask
    call load_original_ortcent !input the particle parameters
    call print_parameters
    if(proc3d%subrec.eq.'y') then
        call reconstruct_substep
    else
        call reconstruct
    end if
    deallocate(particle,ctf2d)
!3D reconstruction    
    


999 end


