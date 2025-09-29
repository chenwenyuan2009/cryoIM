module FFTW3SPACE
      use,intrinsic::iso_c_binding
      use omp_lib
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


  type FFT_plan_cmplx
    integer*8  FFTW_PLAN_cmplx_FWD
    integer*8  FFTW_PLAN_cmplx_BCK
    logical    intialcmplx
    integer nd(3)
  end type


   integer max_plan_number
   parameter(max_plan_number=10)

   type(FFT_plan_cmplx)::FFT_plan(max_plan_number)
   type(FFT_plan_cmplx)::local_FFT_plan

   integer nthread

   integer(4)::fftw_status

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

   interface find_FFT_Plan
      module procedure look_up_plan
   end interface
contains

subroutine new(nthreads)
   integer nthreads,i
   do i=1,max_plan_number
     FFT_plan(i)%FFTW_PLAN_cmplx_BCK=0
     FFT_plan(i)%FFTW_PLAN_cmplx_FWD=0
     FFT_plan(i)%intialcmplx=.FALSE.
     FFT_plan(i)%nd(:)=0
   end do
   nthread=nthreads

   local_FFT_plan%FFTW_PLAN_cmplx_BCK=0
   local_FFT_plan%FFTW_PLAN_cmplx_FWD=0
   local_FFT_plan%intialcmplx=.FALSE.
   local_FFT_plan%nd(:)=0

   call sfftw_init_threads(fftw_status)
end subroutine

subroutine  intialplanCmplx1d(fin,fout,n)
            integer :: n
            complex fin(N),fout(N)
            if(.not.cache_in_plan(n,1,1)) call look_up_plan(fin,fout,n,1,1)
end subroutine  intialplanCmplx1d


subroutine  FFTcmplx1d_f(fin,fout,n)
           implicit none
           integer :: n
           complex fin(N),fout(N)
           if(.not.cache_in_plan(n,1,1)) call look_up_plan(fin,fout,n,1,1)
            call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_FWD, fin,fout)
            fout=fout/float(n)
end subroutine  FFTcmplx1d_f

subroutine  FFTcmplx1d_B(fout,fin,n)
            implicit none
            integer :: n
            complex fin(N),fout(N)
            if(.not.cache_in_plan(n,1,1)) call look_up_plan(fin,fout,n,1,1)
            call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_BCK, fout,fin)
end subroutine  FFTcmplx1d_B



subroutine  intialplanCmplx2d(fin,fout,m,n)
           implicit none
           integer :: m,n
           complex fin(N,M),fout(N,M)
           if(.not.cache_in_plan(m,n,1))  call look_up_plan(fin,fout,m,n,1)
end subroutine  intialplanCmplx2d



subroutine  FFTcmplx2d_f(fin,fout,m,n)
            implicit none
            integer :: m,n
            complex fin(N,M),fout(N,M)
            if(.not.cache_in_plan(m,n,1)) call look_up_plan(fin,fout,m,n,1)
            call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_FWD, fin,fout)
            fout=fout/float(n*m)
end subroutine  FFTcmplx2d_f

subroutine  FFTcmplx2d_B(fout,fin,m,n)
            implicit none
            integer :: m,n
            complex fin(N,M),fout(N,M)
            if(.not.cache_in_plan(m,n,1)) call look_up_plan(fin,fout,m,n,1)
            call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_BCK,fout,fin)
end subroutine  FFTcmplx2d_B


subroutine  intialplanCmplx3d(fin,fout,m,n,k)
            implicit none
            integer :: m,n,k
            complex fin(N,M,k),fout(N,M,k)
            if(.not.cache_in_plan(m,n,k)) call look_up_plan(fin,fout,m,n,k)
end subroutine  intialplanCmplx3d

subroutine  FFTcmplx3d_f(fin,fout,m,n,k)
           implicit none
           integer :: m,n,k
           complex fin(N,M,k),fout(N,M,k)
           if(.not.cache_in_plan(m,n,k)) call look_up_plan(fin,fout,m,n,K)
           call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_FWD, fin,fout)
           fout=fout/float(n*m*k)
end subroutine  FFTcmplx3d_f

subroutine  FFTcmplx3d_B(fout,fin,m,n,k)
            implicit none
            integer :: m,n,k
            complex fin(N,M,k),fout(N,M,k)
            if(.not.cache_in_plan(m,n,k)) call look_up_plan(fin,fout,m,n,K)
            call sfftw_execute_dft(local_FFT_plan%FFTW_PLAN_cmplx_BCK, fout,fin)
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
    real  density(nx,ny)
 !$omp parallel default(shared) private(hx,ky)
 !$omp do schedule(static)
         do ky=1,ny
          do hx=1,nx
           if(mod((hx-1+ky-1),2).eq.0) cycle
           density(hx,ky)= -density(hx,ky)
        end do
       end do
!$omp end do nowait
!$omp end parallel
end subroutine

subroutine shift3dReal(density,nx,ny,nz)
    use omp_lib
    implicit none
    integer hx,ky,lz,nx,ny,nz
    real  density(nx,ny,nz)
 !$omp parallel default(shared) private(hx,ky,lz)
 !$omp do schedule(static)
        do lz=1,nz
          do ky=1,ny
            do hx=1,nx
             if(mod((hx-1+ky-1+lz-1),2).eq.0) cycle
             density(hx,ky,lz)=-density(hx,ky,lz)
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
    Complex  density(nx,ny)
 !$omp parallel default(shared) private(hx,ky)
 !$omp do schedule(static)
        do hx=1,nx
          do ky=1,ny
            if(mod((hx-1+ky-1),2).eq.0) cycle
            density(ky,hx)=density(ky,hx)*(-1.0)
         end do
       end do
!$omp end do nowait
!$omp end parallel
end subroutine

subroutine shift3dComplex(density,nx,ny,nz)
    use omp_lib
    implicit none
    integer hx,ky,lz,nx,ny,nz
    complex  density(nx,ny,nz)
 !$omp parallel default(shared) private(hx,ky,lz)
 !$omp do schedule(static)
      do lz=1,nz
          do ky=1,ny
            do hx=1,nx
             if(mod((hx-1+ky-1+lz-1),2).eq.0) cycle
             density(hx,ky,lz)=-density(hx,ky,lz)
        end do
       end do
       end do
!$omp end do nowait
!$omp end parallel
end   subroutine

subroutine look_up_plan(fin,fout,nx,ny,nz)
   implicit none
   integer nx,ny,nz,pid,findPid,no_use_plan
   complex  fin(nx,ny,nz),fout(nx,ny,nz)
   logical find_no_use

   findPid=0
   no_use_plan=0
   find_no_use=.false.

   do pid=1,max_plan_number
      if(FFT_plan(pid)%intialcmplx) then
        if(FFT_plan(pid)%nd(1).eq.nx .and.  &
         FFT_plan(pid)%nd(2).eq.ny .and.  &
         FFT_plan(pid)%nd(3).eq.nz) then
         findPid=pid
         exit
         end if
      else if(.NOT.find_no_use) then
        no_use_plan=pid
        find_no_use=.TRUE.
      end if
  end do

  if(findPid.ne.0) then
     local_FFT_plan=FFT_plan(findPid)
     return
  end if

  if(no_use_plan.eq.0) then
    no_use_plan=1
    if(FFT_plan(no_use_plan)%intialcmplx) then
      call sfftw_destroy_plan(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_FWD)
      call sfftw_destroy_plan(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_BCK)
      FFT_plan(no_use_plan)%nd(1)=0
      FFT_plan(no_use_plan)%nd(2)=0
      FFT_plan(no_use_plan)%nd(3)=0
      FFT_plan(no_use_plan)%intialcmplx=.FALSE.
    end if
  end if

  call sfftw_plan_with_nthreads(nthread)

  if(ny.eq.1 .and. nz.eq.1) then
      CALL SFFTW_PLAN_DFT_1D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_FWD,NX,fin,fout,FFTW_FORWARD, FFTW_MEASURE)
      CALL SFFTW_PLAN_DFT_1D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_BCK,NX,fout,fin,FFTW_BACKWARD,FFTW_MEASURE)
  else if(nz.eq.1) then
      CALL sFFTW_PLAN_DFT_2D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_FWD,NX,NY,fin,fout,FFTW_FORWARD, FFTW_MEASURE)
      CALL sFFTW_PLAN_DFT_2D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_BCK,NX,NY,fout,fin,FFTW_BACKWARD,FFTW_MEASURE)
  else
      CALL SFFTW_PLAN_DFT_3D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_FWD,nx,ny,nz,fin,fout,FFTW_FORWARD, FFTW_ESTIMATE)
      CALL SFFTW_PLAN_DFT_3D(FFT_plan(no_use_plan)%FFTW_PLAN_cmplx_BCK,nx,ny,nz,fout,fin,FFTW_BACKWARD,FFTW_ESTIMATE)
  end if

  FFT_plan(no_use_plan)%nd(1)=nx
  FFT_plan(no_use_plan)%nd(2)=ny
  FFT_plan(no_use_plan)%nd(3)=nz
  FFT_plan(no_use_plan)%intialcmplx=.TRUE.
  local_FFT_plan=FFT_plan(no_use_plan)
end subroutine


subroutine destoryFFT1d(n)
   implicit none
   integer n,pid
   do pid=1,max_plan_number
     if(FFT_plan(pid)%intialcmplx) then
        if(FFT_plan(pid)%nd(1).eq.n .and. FFT_plan(pid)%nd(2).eq.1 .and.FFT_plan(pid)%nd(3).eq.1) then
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_FWD)
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_BCK)
          FFT_plan(pid)%nd(1)=0
          FFT_plan(pid)%nd(2)=0
          FFT_plan(pid)%nd(3)=0
          FFT_plan(pid)%intialcmplx=.FALSE.
        end if
      end if
   end do

end subroutine

subroutine destoryFFT2d(m,n)
  implicit none
   integer m,n,pid
   do pid=1,max_plan_number
     if(FFT_plan(pid)%intialcmplx) then
        if(FFT_plan(pid)%nd(1).eq.m .and. FFT_plan(pid)%nd(2).eq.n .and.FFT_plan(pid)%nd(3).eq.1) then
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_FWD)
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_BCK)
          FFT_plan(pid)%nd(1)=0
          FFT_plan(pid)%nd(2)=0
          FFT_plan(pid)%nd(3)=0
          FFT_plan(pid)%intialcmplx=.FALSE.
        end if
      end if
   end do
end subroutine

subroutine destoryFFT3d(m,n,k)
   implicit none
   integer m,n,k,pid
   do pid=1,max_plan_number
     if(FFT_plan(pid)%intialcmplx) then
        if(FFT_plan(pid)%nd(1).eq.m .and. FFT_plan(pid)%nd(2).eq.n .and.FFT_plan(pid)%nd(3).eq.k) then
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_FWD)
          call sfftw_destroy_plan(FFT_plan(pid)%FFTW_PLAN_cmplx_BCK)
          FFT_plan(pid)%nd(1)=0
          FFT_plan(pid)%nd(2)=0
          FFT_plan(pid)%nd(3)=0
          FFT_plan(pid)%intialcmplx=.FALSE.
        end if
      end if
   end do
end subroutine

logical function cache_in_plan(m,n,k)
  implicit none
  integer m,n,k
  logical isCache
  isCache=.false.
  if(local_FFT_plan%intialcmplx) then
   if(local_FFT_plan%nd(1).eq.m .and. &
      local_FFT_plan%nd(2).eq.n .and. &
      local_FFT_plan%nd(3).eq.k) then
    isCache=.TRUE.
   end if
  end if
  cache_in_plan=isCache
end function

logical function Plan_in_use(pid,nx,ny,nz)
    implicit none
    logical in_use
    integer nx,ny,nz
    integer pid
    in_use=.FALSE.
    if(pid.le.max_plan_number .and. &
       FFT_plan(pid)%intialcmplx .and.  &
       FFT_plan(pid)%nd(1).eq.nx .and.  &
       FFT_plan(pid)%nd(2).eq.ny .and.  &
       FFT_plan(pid)%nd(3).eq.nz) then
       in_use=.TRUE.
    end if
    Plan_in_use=in_use
end function

End module FFTW3SPACE
