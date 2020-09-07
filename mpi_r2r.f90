Program rupak

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'mpif.h'

  integer (C_INTPTR_T), parameter :: nx = 7
  integer (C_INTPTR_T), parameter :: ny = 5
  integer (C_INTPTR_T), parameter :: nh = ( nx / 2 ) + 1

  include "fftw3-mpi.f03"

  integer(C_INTPTR_T) :: i, j
  real(C_DOUBLE), dimension(nx,ny) :: in
  complex(C_DOUBLE_COMPLEX), dimension(nh,ny) :: out,out1
  integer :: ierr, myid, nproc, seed
  type(C_PTR) :: plan, plan1, cdatar, cdatac
  integer(C_INTPTR_T) :: alloc_local, local_ny, local_j_offset
  real(C_DOUBLE), pointer :: idata(:,:)
  complex(C_DOUBLE_complex), pointer :: odata(:,:)

  seed = 123456789
  
  call mpi_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

  do i = 1, nx
    do j = 1, ny
      in(i,j) = ran(seed)
      seed=in(i,j)*seed
       write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, in(i,j)
    end do
  end do
       WRITE ( *, * ) "_____________"

  call fftw_mpi_init()
  
  ! get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d(ny, nh, MPI_COMM_WORLD, &
                                       local_ny, local_j_offset)
  cdatar = fftw_alloc_real(2 * alloc_local)
  cdatac = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cdatar, idata, [2*nh,local_ny])
  
  call c_f_pointer(cdatac, odata, [nh,local_ny])

  ! create MPI plan for out-of-place DFT (note dimension reversal)
  plan = fftw_mpi_plan_dft_r2c_2d(ny, nx, idata, odata, MPI_COMM_WORLD, &
                                                             FFTW_MEASURE)
  plan1 = fftw_mpi_plan_dft_c2r_2d(ny, nx, odata, idata, MPI_COMM_WORLD, &
                                                             FFTW_MEASURE)
  
  ! initialize data to some function my_function(i,j)
  do i = 1, nx
    do j = 1, local_ny
    idata(i, j) = in(i, j + local_j_offset)
    !write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i,j+local_j_offset,idata(i,j),local_j_offset!,myid
    end do
  end do

  ! compute transform (as many times as desired)
  call fftw_mpi_execute_dft_r2c(plan, idata, odata)

  ! COPY TO OUTPUT ARRAY
  do i = 1, nh
    do j = 1, local_ny
      out(i,j+local_j_offset) = odata(i,j)
      !write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j+local_j_offset, odata(i,j)!,out(i,j+local_j_offset)!, local_j_offset
    end do
  end do
  
  call fftw_mpi_execute_dft_c2r(plan1, odata, idata)

  ! COPY TO OUTPUT ARRAY
  do i = 1, nx
    do j = 1, local_ny
      out1(i,j+local_j_offset) = idata(i,j)
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j+local_j_offset, idata(i,j)/(ny*nx)!,out1(i,j+local_j_offset)!, local_j_offset
    end do
  end do
 
  ! deallocate and destroy plans  
  call fftw_destroy_plan(plan)
  call fftw_mpi_cleanup()
  call fftw_free(cdatar)
  call fftw_free(cdatac)

  do i = 1, nh
    do j = 1, ny
      !write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
    end do
  end do

  do i = 1, nx
    do j = 1, ny
      !write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out1(i,j)/(ny*nx)
    end do
  end do

  call mpi_finalize(ierr)

end
