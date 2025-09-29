
program main
    use common_icos_loc_box
    implicit none
    integer iargc,n,k1,i,j,particleID
    real theta,phi,omega,x0,y0,a,b,c
    character(len=54)::char0


    n=iargc()
    if (n.lt.4) then
      call help
      goto 999
    end if
    do i=1,n
        call getarg(i,par(i))
    end do

    call input_parameters
    do i=1,n
 !       print*,trim(coeff(j))
        k1=index(par(i),'=')
        do j=1,totalpars
            if(par(i)(1:k1-1).eq.trim(coeff(j))) then
                if(trim(coeff(j)).eq.'icos_ort') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) icos_ort
                end if
                if(trim(coeff(j)).eq.'local_ort') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) local_ort
                end if
                if(trim(coeff(j)).eq.'local_imgstck') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) local_imgstck
                    proc2d_local_imgstck='y'
                end if
                if(trim(coeff(j)).eq.'local_x') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) loc_x
                end if
                if(trim(coeff(j)).eq.'local_y') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) loc_y
                end if
                if(trim(coeff(j)).eq.'local_z') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) loc_z
                end if
                if(trim(coeff(j)).eq.'apix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) apix
                end if
                if(trim(coeff(j)).eq.'PR_threshold') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) PR_threshold
                end if
                if(trim(coeff(j)).eq.'sub_FFTsize') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sub_FFTsize
                end if
                if(trim(coeff(j)).eq.'icos_imgstck') then
                    proc2d_icos_imgstck='y'
                    read(par(i)(k1+1:len(trim(par(i)))),*) icos_imgstck
                end if
                if(trim(coeff(j)).eq.'sym') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym 
                end if
                if(trim(coeff(j)).eq.'sym_matrix') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) sym_matrix
                    proc2d_sym_matrix='y'
                end if
                if(trim(coeff(j)).eq.'first') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) first
                end if
                if(trim(coeff(j)).eq.'boundX') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) boundX
                end if
                if(trim(coeff(j)).eq.'last') then
                    read(par(i)(k1+1:len(trim(par(i)))),*) last
                    proc2d_last='y'
                end if
                goto 90
            end if
        end do
        call help
        goto 999
90 end do
   
    if((trim(sym).eq.'I2').or.(trim(sym).eq.'i2')) Nsym=60
    if((trim(sym).eq.'I3').or.(trim(sym).eq.'i3')) Nsym=60
    if((trim(sym).eq.'I5').or.(trim(sym).eq.'i5')) Nsym=60
    if((trim(sym).eq.'c1').or.(trim(sym).eq.'C1')) Nsym=1
    if((trim(sym).eq.'c2').or.(trim(sym).eq.'C2')) Nsym=2
    if((trim(sym).eq.'c3').or.(trim(sym).eq.'C3')) Nsym=3
    if((trim(sym).eq.'c4').or.(trim(sym).eq.'C4')) Nsym=4
    if((trim(sym).eq.'c5').or.(trim(sym).eq.'C5')) Nsym=5
    if((trim(sym).eq.'c6').or.(trim(sym).eq.'C6')) Nsym=6
    if((trim(sym).eq.'c7').or.(trim(sym).eq.'C7')) Nsym=7
    if((trim(sym).eq.'c8').or.(trim(sym).eq.'C8')) Nsym=8
    if((trim(sym).eq.'c9').or.(trim(sym).eq.'C9')) Nsym=9
    if((trim(sym).eq.'c10').or.(trim(sym).eq.'C10')) Nsym=10
    if((trim(sym).eq.'c11').or.(trim(sym).eq.'C11')) Nsym=11
    if((trim(sym).eq.'c12').or.(trim(sym).eq.'C12')) Nsym=12
    if((trim(sym).eq.'c13').or.(trim(sym).eq.'C13')) Nsym=13
    if((trim(sym).eq.'c14').or.(trim(sym).eq.'C14')) Nsym=14
    if((trim(sym).eq.'c15').or.(trim(sym).eq.'C15')) Nsym=15
    
   
   if(proc2d_sym_matrix.eq.'y') then
        Nsym=0
        open(100,file=trim(sym_matrix),status='old')
200     read(100,*,end=300) a,b,c
        read(100,*,end=300) a,b,c
        read(100,*,end=300) a,b,c
        Nsym=Nsym+1
        goto 200
300     close(100)
    end if
    allocate(sym_rot(3,3,Nsym))
    sym_rot=.0
    if(proc2d_sym_matrix.eq.'y') then
        call read_sym_matrixa(sym_matrix,sym_rot,Nsym)
        goto 400
    end if    
    if((trim(sym).eq.'i2').or.(trim(sym).eq.'I2')) then
        call icos_matrix_i2(sym_rot,Nsym)
    else if((trim(sym).eq.'i3').or.(trim(sym).eq.'I3')) then
        call icos_matrix_i3(sym_rot,Nsym)
    else if((trim(sym).eq.'i5').or.(trim(sym).eq.'I5')) then
        call icos_matrix_i5(sym_rot,Nsym)
    else
        call C_matrix(sym_rot,Nsym)
    end if

400 call load_ort
   
    if(proc2d_local_imgstck.eq.'y') call loc_box
    if(proc2d_local_imgstck.eq.'n') call loc_box_coordinate
999 end


