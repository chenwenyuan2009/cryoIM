!This package is used to data process for icosahedral virus.
!five modules:1.3D reconstruction based on direct Fourier transform
!             2.random initial model. 3D reconstruct with random orientation and center (0,0).
!             3.calculate particle center according the 3D model according cross-correlation coefficient.
!             4.globe search particle orientation according to the cross-common lines of icosahedron.
!             5.local refine particle orientation and center according to the cross-common lines of icosahedron.
!Hongrong Liu. 04/01/2016.  HUNNU

program main
    use common_icosproc
    implicit none
    integer iargc,n,k1,i,j,particleID
    real theta,phi,omega,x0,y0


    n=iargc()
    if (n.lt.2) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
    end do
    capsidmodel='capsid.mrc'
    coremodel='core.mrc'
    capsid_innerR=0
    capsid_outerR=270.5
    core_radius=250
    PR_threshold_h=1.0
    loc_box=0
    particle_radius=300

    call input_parameters
    do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'mode') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) mode
                end if
                if(trim(coeff(j)).eq.'model3d') then
                    proc3d%model3d='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) model3d
                end if
                if(trim(coeff(j)).eq.'core3d') then
                    proc3d%core3d='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) core3d
                end if
                if(trim(coeff(j)).eq.'templatestck') then
                    proc3d%templatestck='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) templatestck
                end if
                if(trim(coeff(j)).eq.'ortfile') then
                    proc2d%ortfile='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) ortfile
                end if
                if(trim(coeff(j)).eq.'imgstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgstck
                    proc2d%imgstck='y'
                end if
                if(trim(coeff(j)).eq.'newimgstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) newimgstck
                    proc2d%newimgstck='y'
                end if
                if(trim(coeff(j)).eq.'templateort') then
                    proc3d%templateort='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) templateort
                    totaltemplate=0
                    open(1,file=trim(templateort),status='old')
80                  read(1,*,end=85) particleID,theta,phi,omega,x0,y0
                    totaltemplate=totaltemplate+1
                    goto 80
85                  close(1)
                end if

                goto 90
            end if
        end do
        call help
        goto 999
90 end do


    call initialvalue
    do i=1,n
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'minRes') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) minRes
                end if
                if(trim(coeff(j)).eq.'maxRes') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) maxRes
                end if
                if(trim(coeff(j)).eq.'imgmask') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgmask
                    if(imgmask.gt.0) proc2d%imgmask='y'
                end if
                if(trim(coeff(j)).eq.'applyCTF') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%applyCTF
                end if
                if(trim(coeff(j)).eq.'apix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) apix
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
                if(trim(coeff(j)).eq.'first') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) first
                end if
                if(trim(coeff(j)).eq.'last') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) last
                end if
                if(trim(coeff(j)).eq.'PR_threshold') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) PR_threshold
                end if
                if(trim(coeff(j)).eq.'PR_threshold_h') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) PR_threshold_h
                end if
                if(trim(coeff(j)).eq.'sym') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym
                end if
                if(trim(coeff(j)).eq.'Ncycle') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) Ncycle
                end if
                if(trim(coeff(j)).eq.'result3d') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) result3d
                end if
                if(trim(coeff(j)).eq.'newortfile') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) newortfile
                end if
                if(trim(coeff(j)).eq.'searchstep') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) searchstep
                end if
                if(trim(coeff(j)).eq.'realspaceavg') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%realspaceavg
                end if
                if(trim(coeff(j)).eq.'recISAF') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%recISAF
                end if
                if(trim(coeff(j)).eq.'imask') then
                    proc3d%imask='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) imask
                end if
                if(trim(coeff(j)).eq.'subrec') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%subrec
                end if
                if(trim(coeff(j)).eq.'FSC') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%FSC
                end if
                if(trim(coeff(j)).eq.'FRC') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%FRC
                end if
                if(trim(coeff(j)).eq.'capsid_innerR') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) capsid_innerR
                end if
                if(trim(coeff(j)).eq.'capsid_outerR') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) capsid_outerR
                end if
                if(trim(coeff(j)).eq.'core_radius') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) core_radius
                end if
                if(trim(coeff(j)).eq.'capsidmodel') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) capsidmodel
                end if
                if(trim(coeff(j)).eq.'coremodel') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) coremodel
                end if
                if(trim(coeff(j)).eq.'geneFFTsize') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) geneFFTsize
                    proc2d%geneFFTsize='y'
                end if
                if(trim(coeff(j)).eq.'clipsize')then
                    read(par(i)(k1+1:len(trim(par(i)))),*) clipsize
                    proc2d%clipsize='y'
                endif
                if(trim(coeff(j)).eq.'pr_cc')then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%pr_cc
                endif
                if(trim(coeff(j)).eq.'outputprojstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%outputprojstck
                end if
                if(trim(coeff(j)).eq.'outputrawimgstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%outputrawimgstck
                end if
                if(trim(coeff(j)).eq.'centshift') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%centshift
                end if
                if(trim(coeff(j)).eq.'centbound') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) centbound
                end if
                if(trim(coeff(j)).eq.'FFTsize') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) FFTsize
                end if
                if(trim(coeff(j)).eq.'loc_box') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) loc_box
                end if
                if(trim(coeff(j)).eq.'particle_radius') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) particle_radius
                end if
                if(trim(coeff(j)).eq.'PAD') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) Padsize
                end if

                goto 110
            end if
        end do
        call help
        goto 999
110 end do



    call calc_initialvalues



    if(mode.eq.1) then
        call get_randort_fromicosort
        call reconstruct
    end if

    if(mode.eq.5) then
        call calc_center_imgstck
        call output_ort
    end if

    if(mode.eq.0) call genomeimgstck

    if(mode.eq.2) call reconstruct

    if(mode.eq.3) call algin_core

    if(mode.eq.4) call mrcs2stck

    if(mode.eq.6) call algin_core_phase

    if(mode.eq.7) call algin_vertex

    if(mode.eq.8) call algin_findvertex

    if(mode.eq.9) call algin_core_loc

    if(mode.eq.10) call loc_search

    if(mode.eq.11) call phase_redsiual

    if(mode.eq.12) call mapmask

    if(mode.eq.13) call imagepick

    if(mode.eq.14) call tail_mask

    if(mode.eq.15) call test_applyctf

999 end


