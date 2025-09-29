
    subroutine maphead
    use common_merge_icosprocess_subrec
    implicit none
    FFTsize=FFTsize-8
    mrc%nx=FFTsize  !        /* # of columns ( fastest changing in the map    */
    mrc%ny=FFTsize          !  /* # of rows                                     */
    mrc%nz=FFTsize           ! /* # of sections (slowest changing in the map    */

    mrc%mode=2         ! /* data type
              !    0 = image data in bytes
              !   1 = image data in short integer
              !    2 = image data in floats
              !    3 = complex data in complex short integers
              !    4 = complex data in complex reals          */

    mrc%nxstart=-FFTsize/2      !/* number of first column in map (default = 0)   */
    mrc%nystart=-FFTsize/2       !/* number of first row in map (default = 0)      */
    mrc%nzstart=-FFTsize/2       !/* number of first ssection in map (default = 0) */

    mrc%mx=FFTsize-1           !/* number of intervals along X                   */
    mrc%my=FFTsize-1            !/* number of intervals along Y                   */
    mrc%mz=FFTsize-1            !/* number of intervals along Z                   */

    mrc%a=mrc%mx*apix           !/* cell dimensions in X (angstrom)               */
    mrc%b=mrc%my*apix           !/* cell dimensions in Y (angstrom)               */
    mrc%c=mrc%my*apix           !/* cell dimensions in Z (angstrom)               */

    mrc%alpha=90        !/* cell angles between Y and Z                   */
    mrc%beta=90          !/* cell angles between X and Z                   */
    mrc%gamma=90         !/* cell angles between X and Y                   */

    mrc%mapc=1          !/* number of axis corresponding to columns (X)   */
    mrc%mapr=2          !/* number of axis corresponding to rows (Y)      */
    mrc%maps=3          !/* number of axis corresponding to sections (Z)  */

    mrc%amin=0         ! /* minimum density value                         */
    mrc%amax=0          !/* maximum density value                         */
    mrc%amean=0        ! /* mean density value                            */

    mrc%ispg=0         !/* space group number (0 for images)          */
    mrc% nsymbt=0       ! /* # of bytes for symmetry operators             */
    mrc%extra(:)=0     !/* user defined storage space                    */

    mrc%xorigin=-1.0*mrc%nx       !/* X phase origin                                */
    mrc%yorigin=-1.0*mrc%ny       !/* Y phase origin                                */
    mrc%zorigin=-1.0*mrc%nz       !/* z phase origin                                */

    mrc%map(:)='map '
    mrc%machinestamp(:)='lhr'
    mrc%arms=0

    mrc%nlabl=0         !/* # of labels being used in the MRC header      */

    mrc%label(:,:)='Hongrong Liu' !/* actual text labels

    return
    end



