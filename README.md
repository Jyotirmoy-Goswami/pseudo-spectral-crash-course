# Fluid teaching
a repository for teaching algorithms for fluid codes!

Lecture 1: fftw_1d.f95 + burgulence.f95 (test.f95 + test_1.f95)
Lecture 2: fourier_transform.f95 + burgulence.f95 (fourier.f95 + test_1.f95)
Lecture 3: fftw_2d.f95 (2d_fft.f95)
Lecture 4: poisson_2d.f95 (2d_fft.f95)
Lecture 5: fluid_2d.f95 (2d_fft.f95)
Lecture 6: fluid_2d.f95 (2d_fft.f95) + KH instability
Lecture 7: Revision
Lecture 8: omp_do.f95 + omp_nested_do.f95 (omp.f95)
Lecture 9: mpi_do.f95 (mpi.f95)
Lecture 10: mpi_nested_do.f95 + hybrid_do.f95 (mpi.f95 + hybrid.f95)
Lecture 11: fftw_omp.f95 + fftw_mpi.f95 (2d_fft.f95)
Lecture 12: fluid_2d_with_tracers.f95 (2d_fft.f95)


fftw_1d.f95: This code takes a set of random numbers as input, takes their Fourier transform using one dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
burgulence.f95: This code solves one dimensional Burgers equation using pseudo-spectral method for spatial discretization and Adams-Bashforth algorithm for temporal update, with sin wave as initial condition. A shock appears at time t = 1.
fourier_transform.f95: This code is a one dimensional serial Fourier transform solver. No external library is used in the code. It is for teaching / demonstration purpose only explicitly documenting the formulae of Fourier series from standard texts.
fftw_2d.f95: This code takes a set of random numbers as input of two dimensional array, takes their Fourier transform using two dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
poisson_2d.f95: This code evaluates the Poisson equation in two dimension using psedo-spectral method.
fluid_2d.f95: This code evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. The initial condition is taken as counter-streaming flows thereby akin to Kelvin-Helmholtz type instability.
omp_do.f95: This is a program to teach OpenMP do loops.
omp_nested_do.f95: This is a program to teach OpenMP nested-do loops.
mpi_do.f95: This is a program to teach MPI do loops.
mpi_nested_do.f95: This is a program to teach MPI nested-do loops.
hybrid_do.f95: This is a program to teach Hybrid (= OpenMP + MPI) nested-do loops.
fftw_omp.f95: This is a program to teach how to use FFTW routine with OpenMP parallelization. This is a multi-core extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array in multiple cores, takes their Fourier transform using OpenMP parallel two dimensional FFTW library and then calculates the multi-core inverse Fourier transform, thus finally comparing the result with the input random numbers.
fftw_mpi.f95: This is a program to teach how to use FFTW routine with MPI parallelization. This is a multi-node extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array, takes their Fourier transform using MPI parallel two dimensional FFTW library and then calculates the multi-node inverse Fourier transform, thus finally comparing the result with the input random numbers.
fluid_2d_with_tracers.f95: This code is an extension of the code fluid_2d.f95. Identical to the fluid_2d.f95, it also evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. However, this code has several passive tracer particles sprinkled in the fluid and the evolution of the passive tracer particles are followed using Cloud-In-Cell algorithm.


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
