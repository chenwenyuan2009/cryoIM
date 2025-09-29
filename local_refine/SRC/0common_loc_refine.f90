


module common_loc_refine
    implicit none
    integer totalpar,FFTsize0,FFTsize,first,last,maxFR,Ncycle,totalparticles
    integer avg_norm,dev_norm,halfw,subrang,submult,Nsym
    real pi,minRes,maxRes,apix,imgmask,deg2rad,refine_angle,refine_shift,convF
    parameter(totalpar=22)
    parameter(pi=3.14159265359)
    parameter(deg2rad=pi/180.0)
    parameter(halfw=6)
    character(len=64):: par(50),coeff(totalpar)
    character(len=256):: model3d,imgstck,local_ort0,local_newort,sign
    character(len=1):: RAM,contrast_flip
    character(len=4):: sym
    type particles
        integer particleID,imgID,imgparticleID,local_imgID,local_hx0,local_ky0
        real theta,phi,omega,x,y,PR,df1,df2,astigang
        real local_theta0,local_phi0,local_omega0,local_x0,local_y0,local_z0,local_PR0
        real local_theta, local_phi, local_omega, local_x, local_y, local_z,local_PR
        real local_dtheta,local_dphi,local_domega,local_dx,local_dy,local_dz,local_dPR
        character(len=64)::stckfile
        character(len=1) ::quality
    end type
    type(particles),allocatable::particle(:)

    type CTFparameters
        real vol,Cs,ampwg,Bfactor
    end type
    type(CTFparameters)::CTFpar

    type process2d
        character applyCTF,imgmask,norm,imgstck,centshift,first,last,lp,refine_angle,refine_shift,particle_df,Zhigh_df
    end type
    type(process2d)::proc2d


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

