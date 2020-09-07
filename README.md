# Fluid teaching
a repository for teaching algorithms for fluid codes!

Lecture 1: fft_1d.f95 + burgulence.f95 (test.f95 + test_1.f95)
Lecture 2: fourier_transform.f95 + burgulence.f95 (fourier.f95 + test_1.f95)
Lecture 3: fft_2d.f95 (2d_fft.f95)
Lecture 4: poisson_2d.f95 (2d_fft.f95)
Lecture 5: fluid_2d.f95 (2d_fft.f95)
Lecture 6: fluid_2d.f95 (2d_fft.f95) + KH instability
Lecture 7: Revision
Lecture 8: omp_do.f95 + omp_nested_do.f95 (omp.f95)
Lecture 9: mpi_do.f95 (mpi.f95)
Lecture 10: mpi_nested_do.f95 + hybrid_do.f95 (mpi.f95 + hybrid.f95)
Lecture 11: fftw_omp.f95 + fftw_mpi.f95 (2d_fft.f95)
Lecture 12: incomp_hydro_with_tracers.f95 (2d_fft.f95)



!Serial Do Loop compilation
gfortran <program_name.f95>; ./a.out
!OpenMP Do Loop compilation
gfortran -fopenmp <program_name.f95>; ./a.out
!MPI Do Loop compilation
mpif90 <program_name.f95>; mpirun -np 4 ./a.out
!Hybrid Do Loop compilation
mpif90 -fopenmp <program_name.f95>; mpirun -np 4 ./a.out

!Serial FFTW compilation
gfortran -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3 -lm; ./a.out
!OpenMP FFTW compilation
gfortran -fopenmp -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_omp -lfftw3 -lm; ./a.out
!MPI FFTW compilation
mpif90 -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_mpi -lfftw3 -lm; mpirun -np 4 ./a.out
!Hybrid FFTW compilation
mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun -quiet -np 4 ./a.out
