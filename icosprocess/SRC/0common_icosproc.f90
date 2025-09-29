
module common_icosproc
    implicit none
    integer FFTsize,first,last,totalimg,totalpars,minR,maxR,mode,Ntemplate,loc_box  !
    integer totaltemplate,Nlin,Ncycle,useparticle   !globe search and refinement
    integer inter_box,inter_w,Nsym  !reconstruction
    integer corrCTF_defocuscycle
    real apix,minRes,maxRes,imgmask,imask,boundX  !
    real deltaAngle,deltaTrans,searchstep !globe search and refinement
    real PR_threshold,maxcentshift !reconstruction
    real corrCTF_defocusstep
    parameter(totalpars=32)
    character(len=128)::model3d,ortfile,templateort,newortfile,result3d,templatestck
    character(len=64) ::coeff(totalpars),par(50)
    character(len=54) ::stckfile,tnffile,imgstck
    character(len=4)  ::sym

    type process2d
        character applyCTF,imgmask,norm,FRC,imgstck,corrCTF,check
    end type
    type(process2d)::proc2d

    type process3d
        character realspaceavg,model3d,templatestck,templateort,recISAF,imask,subrec,FSC
    end type
    type(process3d)::proc3d

    type particles
        integer particleID,imgID,imgparticleID
        real theta0,phi0,omega0,x0,y0,PR0,df1,df2,astigang
        real theta, phi, omega, x, y, PR
        real dtheta,dphi,domega,dx,dy,dPR
        character(len=54)::stckfile
    end type
    type(particles),allocatable::particle(:)

    type templates
        integer particleID
        real theta,phi,omega,x,y
    end type
    type(templates),allocatable::template(:)

    type CTFparameters
        real vol,Cs,ampwg,Bfactor
    end type
    type(CTFparameters)::CTFpar

    type recISAFs
      integer  (kind=4)  :: TotalParticle
      integer  (kind=4)  :: InnerRadius, OuterRadius
      integer  (kind=4)  :: Lmax, MulMax, MatrixEvenOrderMax, MatrixOddOrderMax, AngleStep, SampleStep
      real     (kind=4)  :: SelectParticalAllImage, SelectParticalOneImage
      real     (kind=4)  :: theta_m, phi_m
      integer  (kind=4)  :: NumSamplingPointAU, NumISAFsCoeff, NumISAFs
      integer  (kind=4)  :: NumFitCoefficientEven, NumFitCoefficientOdd
    end type
    type(recISAFs) :: ISAFpar

    real,allocatable::template2d(:,:,:),Fcut(:,:),ctf2d(:)

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


