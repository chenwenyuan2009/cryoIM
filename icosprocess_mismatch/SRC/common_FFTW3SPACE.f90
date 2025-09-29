module FFTW3SPACE
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)


! 1D complex
  integer*8  FFTW_PLAN_cmplx1d_f
  integer*8  FFTW_PLAN_cmplx1d_B
  logical    intialcmplx1d
  integer n1d
! 2D complex
  integer*8  FFTW_PLAN_cmplx2d_f
  integer*8  FFTW_PLAN_cmplx2d_B
  logical    intialcmplx2d
  integer n2d(2)

  ! 3D cmplex
  integer*8  FFTW_PLAN_cmplx3d_f
  integer*8  FFTW_PLAN_cmplx3d_B
  logical   intialcmplx3d
  integer n3d(3)

  interface FFTNEW
    module procedure new
  end interface

  interface intialplan
    module procedure intialplanCmplx1d,intialplanCmplx2d,intialplanCmplx3d
  end interface

  interface doFFT
    module procedure FFTcmplx1d_f,FFTcmplx2d_f,FFTcmplx3d_f
  end interface

   interface doIFFT
    module procedure FFTcmplx1d_B,FFTcmplx2d_B,FFTcmplx3d_B
  end interface

   interface shift
     module procedure shift1dReal,shift2dReal,shift3dReal
     module procedure shift1dCmplex,shift2dComplex,shift3dComplex
   end interface

   interface destoryplan
     module procedure destoryFFT1d,destoryFFT2d,destoryFFT3d
   end interface
contains

subroutine new()
   FFTW_PLAN_cmplx1d_f=0
   FFTW_PLAN_cmplx1d_B=0
   intialcmplx1d      =.False.
   n1d=0
! 2D complex
   FFTW_PLAN_cmplx2d_f=0
   FFTW_PLAN_cmplx2d_B=0
   intialcmplx2d      =.False.
   n2d=0
  ! 3D cmplex
   FFTW_PLAN_cmplx3d_f=0
   FFTW_PLAN_cmplx3d_B=0
   intialcmplx3d      =.False.
   n3d=0
end subroutine

subroutine  intialplanCmplx1d(fin,fout,n)
            integer :: n
            complex fin(N),fout(N)
           if( intialcmplx1d .eqv. .False. ) then
                CALL SFFTW_PLAN_DFT_1D(FFTW_PLAN_cmplx1d_f,N,fin,fout,FFTW_FORWARD, FFTW_MEASURE)
                CALL SFFTW_PLAN_DFT_1D(FFTW_PLAN_cmplx1d_B,N,fout,fin,FFTW_BACKWARD,FFTW_MEASURE)
                intialcmplx1d= .True.
                n1d=n
          end if
end subroutine  intialplanCmplx1d

subroutine  FFTcmplx1d_f(fin,fout,n)
            integer :: n
            complex fin(N),fout(N)
            if(n1d .ne. n .and. intialcmplx1d .eqv. .True. ) call destoryFFT1d(n1)
            call intialplanCmplx1d(fin,fout,n)
            call sfftw_execute_dft(FFTW_PLAN_cmplx1d_f, fin,fout)
            fout=fout/float(n)
end subroutine  FFTcmplx1d_f

subroutine  FFTcmplx1d_B(fout,fin,n)
             integer :: n
            complex fin(N),fout(N)
            if(n1d.ne.n.and.intialcmplx1d .eqv..True. ) call destoryFFT1d(n1)
            call intialplanCmplx1d(fin,fout,n)
            call sfftw_execute_dft(FFTW_PLAN_cmplx1d_B, fout,fin)
end subroutine  FFTcmplx1d_B

subroutine  destoryFFT1d(n1)
         integer n1
         call sfftw_destroy_plan(FFTW_PLAN_cmplx1d_f)
         call sfftw_destroy_plan(FFTW_PLAN_cmplx1d_B)
         intialcmplx1d=.False.
         n1=n1
         n1d=0
end subroutine  destoryFFT1d

subroutine  destoryFFT2d(nx,ny)
         integer nx,ny
         call sfftw_destroy_plan(FFTW_PLAN_cmplx2d_f)
         call sfftw_destroy_plan(FFTW_PLAN_cmplx2d_B)
         intialcmplx2d=.False.
         nx=nx
         ny=ny
         n2d=0
end subroutine  destoryFFT2d

subroutine  destoryFFT3d(nx,ny,nz)
         integer nx,ny,nz
         call sfftw_destroy_plan(FFTW_PLAN_cmplx3d_f)
         call sfftw_destroy_plan(FFTW_PLAN_cmplx3d_B)
         intialcmplx3d=.False.
         nx=nx
         ny=ny
         nz=nz
         n3d=0
end subroutine  destoryFFT3d


subroutine  intialplanCmplx2d(fin,fout,m,n)
           integer :: m,n
           complex fin(N,M),fout(N,M)
           if( intialcmplx2d .eqv. .False.) then
              CALL SFFTW_PLAN_DFT_2D(FFTW_PLAN_cmplx2d_f,N,M,fin,fout,FFTW_FORWARD, FFTW_MEASURE)
              CALL SFFTW_PLAN_DFT_2D(FFTW_PLAN_cmplx2d_B,N,M,fout,fin,FFTW_BACKWARD,FFTW_MEASURE)
              intialcmplx2d=.True.
              n2d(1)=N
              n2d(2)=M
          end if
end subroutine  intialplanCmplx2d

subroutine  FFTcmplx2d_f(fin,fout,m,n)
            integer :: m,n
            complex fin(N,M),fout(N,M)
            if(n2d(1).ne.n.and.n2d(2).ne.m .and.intialcmplx2d .eqv..True. ) then
                call destoryFFT2d(m,n)
            end if
            call intialplanCmplx2d(fin,fout,m,n)
            call sfftw_execute_dft(FFTW_PLAN_cmplx2d_f, fin,fout)
             fout=fout/float(n*m)
end subroutine  FFTcmplx2d_f

subroutine  FFTcmplx2d_B(fout,fin,m,n)
            integer :: m,n
            complex fin(N,M),fout(N,M)
            if(n2d(1).ne.n.and.n2d(2).ne.m .and.intialcmplx2d .eqv..True. ) then
                call destoryFFT2d(m,n)
            end if
            call intialplanCmplx2d(fin,fout,m,n)
            call sfftw_execute_dft(FFTW_PLAN_cmplx2d_B, fout,fin)
end subroutine  FFTcmplx2d_B


subroutine  intialplanCmplx3d(fin,fout,m,n,k)
            integer :: m,n,k
            complex fin(N,M,k),fout(N,M,k)
            if(intialcmplx3d .eqv. .False.) then
                CALL SFFTW_PLAN_DFT_3D(FFTW_PLAN_cmplx3d_f,m,n,k,fin,fout,FFTW_FORWARD, FFTW_ESTIMATE)
                CALL SFFTW_PLAN_DFT_3D(FFTW_PLAN_cmplx3d_B,m,n,k,fout,fin,FFTW_BACKWARD,FFTW_ESTIMATE)
                intialcmplx3d=.True.
                n3d(1)=N
                n3d(2)=M
                n3d(3)=k
           end if
end subroutine  intialplanCmplx3d

subroutine  FFTcmplx3d_f(fin,fout,m,n,k)
           integer :: m,n,k
           complex fin(N,M,k),fout(N,M,k)
           if(n3d(1).ne.n.and.n3d(2).ne.m.and.n3d(3).ne.k   &
             .and.intialcmplx3d .eqv..True. ) then
                call destoryFFT3d(m,n,k)
            end if
           call intialplanCmplx3d(fin,fout,m,n,k)
           call sfftw_execute_dft(FFTW_PLAN_cmplx3d_f, fin,fout)
           fout=fout/float(n*m*k)
end subroutine  FFTcmplx3d_f

subroutine  FFTcmplx3d_B(fout,fin,m,n,k)
            integer :: m,n,k
            complex fin(N,M,k),fout(N,M,k)
            if(n3d(1).ne.n.and.n3d(2).ne.m.and.n3d(3).ne.k   &
             .and.intialcmplx3d .eqv..True. ) then
                call destoryFFT3d(m,n,k)
            end if
            call intialplanCmplx3d(fin,fout,m,n,k)
           call sfftw_execute_dft(FFTW_PLAN_cmplx3d_B, fout,fin)
end subroutine  FFTcmplx3d_B

subroutine shift1dReal(density,n)
    use omp_lib
    implicit none
    integer hx,half,n
    real  density(-n/2:n/2-1)
    half=n/2
 !$omp parallel default(shared) private(hx)
 !$omp do schedule(static)
        do hx=-half,half-1
          density(hx)=density(hx)*(-1.0)**(hx)
        end do
!$omp end do nowait
!$omp end parallel
end subroutine shift1dReal

subroutine shift1dCmplex(FFT,n)
    use omp_lib
    implicit none
    integer hx,half,n
    complex FFT(-n/2:n/2-1)
    half=n/2
 !$omp parallel default(shared) private(hx)
 !$omp do schedule(static)
        do hx=-half,half-1
             FFT(hx)=FFT(hx)*(-1.0)**(hx)
         end do
!$omp end do nowait
!$omp end parallel
end subroutine  shift1dCmplex

subroutine shift2dReal(density,nx,ny)
    use omp_lib
    implicit none
    integer hx,ky,nx,ny
    real  density(-nx/2:nx/2-1,-ny/2:ny/2-1)
 !$omp parallel default(shared) private(hx,ky)
 !$omp do schedule(static)
        do hx=-nx/2,nx/2-1
          do ky=-ny/2,ny/2-1
           density(hx,ky)= density(hx,ky)*(-1.0)**(hx+ky)
        end do
       end do
!$omp end do nowait
!$omp end parallel
end subroutine

subroutine shift3dReal(density,nx,ny,nz)
    use omp_lib
    implicit none
    integer hx,ky,lz,nx,ny,nz
    real  density(-nx/2:nx/2-1,-ny/2:ny/2-1,-nz/2:nz/2-1)
 !$omp parallel default(shared) private(hx,ky,lz)
 !$omp do schedule(static)
      do hx=-nx/2,nx/2-1
          do ky=-ny/2,ny/2-1
            do lz=-nz/2,nz/2-1
             density(hx,ky,lz)=density(hx,ky,lz)*(-1.0)**(hx+ky+lz)
        end do
       end do
       end do
!$omp end do nowait
!$omp end parallel
end   subroutine

subroutine shift2dComplex(density,nx,ny)
    use omp_lib
    implicit none
    integer hx,ky,nx,ny
    Complex  density(-nx/2:nx/2-1,-ny/2:ny/2-1)
 !$omp parallel default(shared) private(hx,ky)
 !$omp do schedule(static)
         do hx=-nx/2,nx/2-1
          do ky=-ny/2,ny/2-1
            density(hx,ky)=density(hx,ky)*(-1.0)**(hx+ky)
         end do
       end do
!$omp end do nowait
!$omp end parallel
end subroutine

subroutine shift3dComplex(density,nx,ny,nz)
    use omp_lib
    implicit none
    integer hx,ky,lz,nx,ny,nz
    complex  density(-nx/2:nx/2-1,-ny/2:ny/2-1,-nz/2:nz/2-1)
 !$omp parallel default(shared) private(hx,ky,lz)
 !$omp do schedule(static)
        do hx=-nx/2,nx/2-1
          do ky=-ny/2,ny/2-1
            do lz=-nz/2,nz/2-1
             density(hx,ky,lz)= density(hx,ky,lz)*(-1.0)**(hx+ky+lz)
        end do
       end do
       end do
!$omp end do nowait
!$omp end parallel
end   subroutine

End module FFTW3SPACE
