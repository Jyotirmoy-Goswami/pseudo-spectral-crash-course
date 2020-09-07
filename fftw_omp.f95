  program fftw_omp
  use omp_lib
  implicit none
  integer ( kind = 4 ), parameter :: nx = 8
  integer ( kind = 4 ), parameter :: ny = 10
  integer ( kind = 4 ), parameter :: nh = ( nx / 2 ) + 1
  include "fftw3.f"
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) thread_num,num_thread,proc_num,iret
  real ( kind = 8 ) in(nx,ny)
  real ( kind = 8 ) in2(nx,ny)
  integer ( kind = 4 ) j
  complex ( kind = 8 ) out(nh,ny)
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 4 ) seed
  seed = 123456789
  proc_num = omp_get_num_procs()
  thread_num = 15
  call omp_set_num_threads (thread_num)
  print*, thread_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate FFTW3 on a 2D real array'
  write ( *, '(a,i8)' ) '  NX = ', nx
  write ( *, '(a,i8)' ) '  NY = ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
  write ( *, '(a)' ) '  the data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'
!
!  Compute the data.
!
  do j = 1, ny
    do i = 1, nx
      in(i,j) = ran(seed)
      seed=in(i,j)*seed
    end do
  end do
write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '
  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, j, in(i,j)
    end do
  end do
!
!  Make a plan for the FFT, and forward transform the data.
!
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(15)
  call dfftw_plan_dft_r2c_2d_ ( plan_forward, nx, ny, in, out, FFTW_ESTIMATE )
  call dfftw_execute_ ( plan_forward )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '
  do i = 1, nh
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
    end do
  end do
!
!  Make a plan for the backward FFT, and recover the original data.
!
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(15)
  call dfftw_plan_dft_c2r_2d_ ( plan_backward, nx, ny, out, in2, FFTW_ESTIMATE )
  call dfftw_execute_ ( plan_backward )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
  write ( *, '(a)' ) ' '
  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) &
        i, j, in2(i,j) / real ( nx * ny, kind = 8 )
    end do
  end do
!
!  Discard the information associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )
!  return
end program fftw_omp
