
module common_proc_icos
     implicit none
    integer FFTsize,sub_FFTsize,first,last,totalpars,totalparticle,Nsym
    real apix,PR_threshold,loc_x,loc_y,loc_z,boundX
    real,allocatable::sym_rot(:,:,:)
    parameter(totalpars=9)
    character(len=64)::ort0,newort,local_imgstck,icos_imgstck,sym_matrix,calc_matrix_3f
    character(len=64)::coeff(totalpars),par(50),calc_randort_2f,calc_matrix_2f
    character(len=64)::calc_matrix_5f,transform_ort_2fTo5f,transform_ort_2fTo3f
    character(len=1)::proc2d_icos_imgstck,proc2d_sym_matrix,proc2d_last,proc2d_local_imgstck
    character(len=2)::sym

    type particles
        integer particleID,imgID,imgparticleID
        real theta0,phi0,omega0,x0,y0,PR0,df1,df2,astigang
        real theta, phi, omega, x, y, PR
        real dtheta,dphi,domega,dx,dy,dPR
        character(len=64)::stckfile
        character(len=1) ::status
    end type
    type(particles),allocatable::particle(:)


    type mrchead
        integer nx  !        /* # of columns ( fastest changing in the map    */
        integer ny          !  /* # of rows                                     */
        integer nz           ! /* # of sections (slowest changing in the map    */

        integer mode         ! /* data type
              !    0 = image data in bytes
              !   1 = image data in short integer
              !    2 = image data in floats
              !    3 = complex data in complex short integers
              !    4 = complex data in complex reals          */

         integer nxstart      !/* number of first column in map (default = 0)   */
         integer nystart       !/* number of first row in map (default = 0)      */
         integer nzstart       !/* number of first ssection in map (default = 0) */

         integer mx           !/* number of intervals along X                   */
         integer my            !/* number of intervals along Y                   */
         integer mz            !/* number of intervals along Z                   */

         real a           !/* cell dimensions in X (angstrom)               */
         real b           !/* cell dimensions in Y (angstrom)               */
         real c           !/* cell dimensions in Z (angstrom)               */

         real alpha        !/* cell angles between Y and Z                   */
         real beta          !/* cell angles between X and Z                   */
         real gamma         !/* cell angles between X and Y                   */

         integer mapc          !/* number of axis corresponding to columns (X)   */
         integer mapr          !/* number of axis corresponding to rows (Y)      */
         integer maps          !/* number of axis corresponding to sections (Z)  */

         real amin         ! /* minimum density value                         */
         real amax          !/* maximum density value                         */
         real amean        ! /* mean density value                            */

         integer ispg         !/* space group number (0 for images)          */
         integer nsymbt       ! /* # of bytes for symmetry operators             */
         integer extra(25)     !/* user defined storage space                    */

         real xorigin       !/* X phase origin                                */
         real yorigin       !/* Y phase origin                                */
         real zorigin       !/* z phase origin                                */

         character  map(4)
         character  machinestamp(4)
         real arms

         integer nlabl         !/* # of labels being used in the MRC header      */

         character label(10,80) !/* actual text labels
    end type
    type(mrchead):: mrc

end module


