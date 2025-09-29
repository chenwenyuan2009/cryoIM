subroutine calc_initialvalues
use common_icosproc
    implicit none

    minRes=apix*FFTsize/minRes
    maxRes=apix*FFTsize/maxRes


    minR  =int(minRes)
    maxR  =int(maxRes)
    if(loc_box.ne.0)then
        Padsize=Padsize*loc_box
    else
        Padsize=Padsize*FFTsize
    endif
!    if(mod(maxR,2).ne.0) maxR=maxR+1
    deltaAngle=360.0/(maxRes*sqrt(5.0)*3.14159); !angle step for ort refine
    deltaTrans=float(FFTsize)/(maxRes*sqrt(5.0)); !cent shift step for cent refine
    if(trim(sym).eq.'c1') then
        Nsym=1
        Nlin=1
    end if
    if(trim(sym).eq.'icos') then
        Nsym=60
        Nlin=60
    end if
    if(trim(sym).eq.'c12') then
        Nsym=12
        Nlin=12
    end if
    if(trim(sym).eq.'c15') then
        Nsym=5
        Nlin=5
    end if
    if(trim(sym).eq.'c6') then
        Nsym=6
        Nlin=6
    end if
    if(proc2d%imgmask.eq.'y') then
        imgmask=imgmask/apix
    else
        imgmask=FFTsize/2-2
    end if

    capsid_innerR=capsid_innerR/apix
    capsid_outerR=capsid_outerR/apix
    core_radius  =core_radius/apix

    allocate(particle(first:last))
    allocate(template(0:totaltemplate))
    allocate(template2d(FFTsize**2,2,0:totaltemplate))
    allocate(Fcut(minR:maxR,0:totaltemplate))
    allocate(ctf2d(FFTsize**2))
    ctf2d=1.0

    if(mode.ne.12 .and. mode.ne.14) call load_original_ortcent !input the particle parameters
    if(proc2d%imgstck.eq.'n') call checkstckfile

    if(mode.eq.1) then
        print*,'random 3D reconstruction:'
        print*,'ortfile =',trim(ortfile)
        print*,'result  =',trim(result3d)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'PR_threshold=',PR_threshold
        print*,'totalparticle:',last-first+1
        print*,'Useparticle=',useparticle
        print*,'sym     =',trim(sym)
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.2) then
        print*,'3D reconstruction parameters:'
        print*,'result  =',trim(result3d)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'PR_threshold=',PR_threshold
        print*,'totalparticle:',last-first+1
        print*,'Useparticle=',useparticle
        print*,'sym     =',trim(sym)
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.3) then
        print*,'calculate particle center according initial model:'
        print*,'ortcent =',trim(ortfile)
        if(proc3d%model3d.eq.'y') then
            print*,'model3d =',trim(model3d)
        else
            print*,'calculate particle center according to rotate 180 degree correlation.'
        end if
        print*,'newcent =',trim(newortfile)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.4) then
        print*,'globe search:'
        print*,'ortcent =',trim(ortfile)
        print*,'model3d =',trim(model3d)
        print*,'temport =',trim(templateort)
        print*,'newortfiel=',trim(newortfile)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'minR    =',minR
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'sym     =',trim(sym)
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.5) then
        print*,'local refinement:'
        print*,'ortcent =',trim(ortfile)
        print*,'model3d =',trim(model3d)
        print*,'temport =',trim(templateort)
        print*,'newortfile=',trim(newortfile)
        print*,'Ncycle  =',Ncycle
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'minR    =',minR
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'sym     =',trim(sym)
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.0) then
        print*,'genomeimgstck subtract:'
        print*,'ortfile =',trim(ortfile)
        print*,'result  =',trim(result3d)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if
        if(mode.eq.7) then
        print*,'vertex refinement:'
        print*,'ortcent =',trim(ortfile)
        print*,'model3d =',trim(model3d)
        print*,'temport =',trim(templateort)
        print*,'newortfile=',trim(newortfile)
        print*,'FFTsize =',FFTsize
        print*,'apix    =',apix
        print*,'minR    =',minR
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'Padsiez=',Padsize
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if



    return
end

