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
    use common_icosproc
    implicit none
    integer iargc,n,k1,i,j,particleID1 !,char0
    real theta1,phi1,omega1,x01,y01
    !character(len=54)::


    n=iargc()
    if (n.lt.2) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
        !print*,par(i)
    end do

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
                if(trim(coeff(j)).eq.'templatestck') then
                    proc3d%templatestck='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) templatestck
                end if
                if(trim(coeff(j)).eq.'ortfile') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ortfile
                end if
                if(trim(coeff(j)).eq.'templateort') then
                    proc3d%templateort='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) templateort
                     totaltemplate=0
                    open(1,file=trim(templateort),status='old')
80                  read(1,*,end=85) particleID1,theta1,phi1,omega1,x01,y01
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

!______________________________________________________________________
! 2d projection

    if(mode.eq.0) then
        if(proc3d%model3d.ne.'y') then
            print*,'should input model3d.'
            stop
        end if
        if(proc3d%templateort.ne.'y') then
            print*,'should input templateort.'
            stop
        end if
        if(proc3d%templatestck.ne.'y') then
            print*,'should input templatestck.'
            stop
        end if
        allocate(template(totaltemplate))
        call project2d
        stop
    end if
!_______________________________________________________________________

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
                if(trim(coeff(j)).eq.'check') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc2d%check
                end if
                if(trim(coeff(j)).eq.'apix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) apix
                end if
                if(trim(coeff(j)).eq.'vol') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ctfpar%vol
                end if
                if(trim(coeff(j)).eq.'Cs') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ctfpar%Cs
                end if
                if(trim(coeff(j)).eq.'ampwg') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ctfpar%ampwg
                end if
                if(trim(coeff(j)).eq.'bfactor') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) ctfpar%bfactor
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
                if(trim(coeff(j)).eq.'recisaf') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%recisaf
                end if
                if(trim(coeff(j)).eq.'subrec') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%subrec
                end if
                if(trim(coeff(j)).eq.'fsc') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) proc3d%fsc
                end if
                if(trim(coeff(j)).eq.'imgstck') then
                    proc2d%imgstck='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) imgstck
                end if
                if(trim(coeff(j)).eq.'sym') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym
                end if
                if(trim(coeff(j)).eq.'maxcentshift') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) maxcentshift
                end if
                if(trim(coeff(j)).eq.'boundX') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) boundX
                end if
                if(trim(coeff(j)).eq.'corrctf') then
                    proc2d%corrctf='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) corrctf_defocusstep,corrctf_defocuscycle 
                end if
                if(trim(coeff(j)).eq.'loc_box') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) loc_box
                end if

                goto 110
            end if
        end do
        call help
        goto 999
110 end do

    call calc_initialvalues

!3d reconstruction    
    
    if((mode.eq.1).and.(proc2d%corrctf.eq.'y')) then
        if((proc3d%subrec.eq.'Y').or.(proc3d%subrec.eq.'y')) then
            call reconstruct_corr_ctf_substep
        else
            call reconstruct_corrctf
        end if
        goto 999
    end if
    if(mode.eq.1) then
        if(trim(sym).ne.'icos') then
            if((proc3d%subrec.eq.'Y').or.(proc3d%subrec.eq.'y')) then
                print*,'perform reconstruct_c1_substep'
                call reconstruct_c1_substep
            else
                call reconstruct_c1
            end if
        else if((proc3d%subrec.eq.'Y').or.(proc3d%subrec.eq.'y')) then
            if((proc3d%fsc.eq.'Y').or.(proc3d%fsc.eq.'y')) then
                tnffile=trim(result3d)//'oddtnf'
                call reconstruct_half
                first=first+1
                tnffile=trim(result3d)//'evetnf'
                call reconstruct_half
            else
                call reconstruct_substep
            end if
        else if((proc3d%fsc.eq.'Y').or.(proc3d%fsc.eq.'y')) then
            tnffile='oddtnf'
            call reconstruct_half
            first=first+1
            tnffile='evetnf'
            call reconstruct_half
            call reconstruct_fsc
        else
            call reconstruct
        end if
    end if

    !random reconstruction
    if(mode.eq.2) then 
        particle(:)%theta=.0;particle(:)%phi=.0;particle(:)%omega=.0;particle(:)%x=.0;particle(:)%y=.0;particle(:)%pr=.0               
        call calc_center_selfrot180
        call get_randomort
        call reconstruct
    end if

    !calculate center according model3d
    if(mode.eq.3) then
        if(proc3d%model3d.eq.'y') call calc_center    
        if(proc3d%model3d.ne.'y') call calc_center_selfrot180
        call output_ort
    end if

    if(mode.eq.4) then
        call load_template
        call ort_globesearch
    end if

    if(mode.eq.5) then
        call load_template
        call ort_refinement
    end if

    if(mode.eq.6) then
        if((proc3d%subrec.eq.'y').or.(proc3d%subrec.eq.'y')) then
        call reconstruct_local_substep
        else
        call reconstruct_local
        endif
    end if
     deallocate(particle,template,template2d,Fcut,ctf2d)
999 end


