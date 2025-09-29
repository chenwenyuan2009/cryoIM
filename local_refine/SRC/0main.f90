
program main
use common_loc_refine
    implicit none
    integer iargc,n
    integer i,j,k1

    n=iargc()
    if (n.lt.5) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
    end do

    call input_parameters

    do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpar
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'local_ort0') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) local_ort0
                end if
                if(trim(coeff(j)).eq.'local_newort') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) local_newort
                end if
                if(trim(coeff(j)).eq.'model3d') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) model3d
                end if
                if(trim(coeff(j)).eq.'imgstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgstck
                    proc2d%imgstck='y'
                end if
                if(trim(coeff(j)).eq.'maxRes') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) maxRes
                end if
                if(trim(coeff(j)).eq.'minRes') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) minRes
                end if
                if(trim(coeff(j)).eq.'apix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) apix
                end if
                if(trim(coeff(j)).eq.'first') then
                    proc2d%first='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) first
                end if
                if(trim(coeff(j)).eq.'last') then
                    proc2d%last='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) last
                end if
                if(trim(coeff(j)).eq.'imgmask') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgmask
                    proc2d%imgmask='y'
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
                if(trim(coeff(j)).eq.'refine_angle') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) refine_angle
                    proc2d%refine_angle='y'
                end if
                if(trim(coeff(j)).eq.'refine_shift') then
                    proc2d%centshift='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) refine_shift
                    proc2d%refine_shift='y'
                end if
                if(trim(coeff(j)).eq.'Ncycle') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) Ncycle
                end if
                if(trim(coeff(j)).eq.'convF') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) convF
                end if
                if(trim(coeff(j)).eq.'particle_df') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%particle_df
                end if
                if(trim(coeff(j)).eq.'Zhigh_df') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%Zhigh_df
                end if
                if(trim(coeff(j)).eq.'sym') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym
                end if
                if(trim(coeff(j)).eq.'sign') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sign
                end if
                goto 110
            end if
        end do
        call help
        goto 999
110 end do

    call calc_initialvilues
    print*,'FFTsize0,FFTsize,maxRes,Ncycle,convF,refine_angle,refine_shift:'
    print*,FFTsize0,FFTsize,maxRes,Ncycle,convF,refine_angle,refine_shift
    if((proc2d%Zhigh_df.eq.'y').or.(proc2d%Zhigh_df.eq.'Y')) then
        call refine_PR_Zhigh
    else if((proc2d%particle_df.eq.'y').or.(proc2d%particle_df.eq.'Y'))then 
        call refine_particle_df
    else
        call refine_PR
    end if
999 end


