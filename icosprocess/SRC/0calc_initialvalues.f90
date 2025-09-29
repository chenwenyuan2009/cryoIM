subroutine calc_initialvalues
use common_icosproc
    implicit none
    if(mode.ne.6)then
    !print*,maxRes
    minRes=apix*FFTsize/minRes
    maxRes=apix*FFTsize/maxRes
    !print*,apix,FFTsize
    else
    minRes=apix*loc_box/minRes
    maxRes=apix*loc_box/maxRes
    endif
    minR  =int(minRes)
    maxR  =int(maxRes)
    

    if(mod(maxR,2).ne.0) maxR=maxR+1
    if(proc2d%imgmask.eq.'y')  imgmask=imgmask/apix

    deltaAngle=360.0/(maxRes*sqrt(5.0)*3.14159); !angle step for ort refine
    deltaTrans=float(FFTsize)/(maxRes*sqrt(5.0)); !cent shift step for cent refine
    if(trim(sym).eq.'c1') then
        Nsym=1
        Nlin=1
    end if
    if(trim(sym).eq.'c2') then
        Nsym=2
        Nlin=2
    end if
    if(trim(sym).eq.'c7') then
        Nsym=7
        Nlin=7
    end if
    if(trim(sym).eq.'c8') then
        Nsym=8
        Nlin=8
    end if
    if(trim(sym).eq.'c12') then
        Nsym=12
        Nlin=12
    end if
    if(trim(sym).eq.'c15') then
        Nsym=15
        Nlin=15
    end if
    if(trim(sym).eq.'c5a') then
        Nsym=5
        Nlin=5
    end if
    if(trim(sym).eq.'c6') then
        Nsym=6
        Nlin=6
    end if
    if(trim(sym).eq.'c3') then
        Nsym=3
        Nlin=3
    end if
    if(trim(sym).eq.'c4') then
        Nsym=4
        Nlin=4
    end if
    allocate(particle(first:last))
    allocate(template(0:totaltemplate))
    allocate(template2d(FFTsize**2,2,0:totaltemplate))
    allocate(Fcut(minR:maxR,0:totaltemplate))
    if(mode.ne.6)then
    allocate(ctf2d(FFTsize**2))
    else
    allocate(ctf2d(loc_box**2))
    endif
    ctf2d=1.0

    call load_original_ortcent !input the particle parameters
    call checkstckfile

    if(mode.eq.1) then
        print*,'3D reconstruction parameters:'
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
        print*,'applyCTF=',proc2d%applyCTF
        print*,'sym=', sym
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    if(mode.eq.2) then
        print*,'random 3D reconstruction:'
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
        print*,'applyCTF=',proc2d%applyCTF
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

        if(mode.eq.6) then
        print*,'local reconstrction:'
        print*,'ortcent =',trim(ortfile)
        print*,'model3d =',trim(model3d)
        print*,'FFTsize =',FFTsize
        print*,'loc_box=',loc_box
        print*,'apix    =',apix
        print*,'minR    =',minR
        print*,'maxR    =',maxR
        print*,'imgmask =',imgmask
        print*,'first   =',first
        print*,'last    =',last
        print*,'applyCTF=',proc2d%applyCTF
        print*,'Useparticle=',useparticle
        if((proc2d%applyCTF.eq.'y').or.(proc2d%applyCTF.eq.'Y')) then
            print*,'vol     =',CTFpar%vol
            print*,'Cs      =',CTFpar%Cs
            print*,'Bfactor =',CTFpar%Bfactor
            print*,'ampwg   =',CTFpar%ampwg
        end if
    end if

    return
end

