Solving PDEs by Pseudo Spectral Method
===============================================
A repository of teaching numerical algorithms for solving fluid equations. The lectures slides can be found [here](PseudoSpectralTeachingSlides.pdf).

Instructor: | [Dr. Rupak Mukherjee](https://github.com/RupakMukherjee)

Contributor: | [Dr. Sayan Adhikari](https://github.com/sayanadhikari)


We have described pseudo-spectral method for spatial discretization. The equations solved in this lecture series are two dimensional incompressible fluid equations in vorticity (<img src="https://render.githubusercontent.com/render/math?math=\omega">) - stream function (<img src="https://render.githubusercontent.com/render/math?math=\psi">) formalism.


<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial \vec{\omega}}{\partial t} %2B \vec{u} \cdot \vec{\nabla} \vec{\omega} = \nu \nabla^2 \vec{\omega}">
<img src="https://render.githubusercontent.com/render/math?math=\vec{\omega} = \vec{\nabla} \times \vec{u}">
<img src="https://render.githubusercontent.com/render/math?math=\nabla^2 \psi = - \omega ">


The wave-form is assumed to be <img src="https://render.githubusercontent.com/render/math?math=\exp(i(kx-wt))">. Hence derivative w.r.t. x gives ik only (and NOT -ik).

A series of lectures (in [Hinglish](https://en.wikipedia.org/wiki/Hinglish)) are available on youtube based on the materials archived in this repository. Click on the link provided below to access the lectures.


[Webinar on "Pseudo Spectral Method"](https://www.youtube.com/playlist?list=PLbX_ZyxeXxSJXnIAnkhhAsIAV-Ld0Awsu)

[![Webinar on "Pseudo Spectral Method"](http://img.youtube.com/vi/m_dle8vr3dU/0.jpg)](https://www.youtube.com/embed/videoseries?list=PLbX_ZyxeXxSJXnIAnkhhAsIAV-Ld0Awsu)

# Lecture description

[Lecture 1](lecture1.md) | [Lecture 2](lecture2.md) | [Lecture 3](lecture3.md) | [Lecture 4](lecture4.md) | [Lecture 5](lecture5.md) | [Lecture 6](lecture6.md)
------------ | ------------- | ------------ | ------------- | ------------ | ------------- 
[Resources](resources/resource_lecture1.md) | [Resources](resources/resource_lecture2.md) | None      | None      | None      | [Resources](resources/resource_lecture6.md) 
None      | None      | [Quiz](quiz/quiz_lecture3.md)      | None      | [Quiz](quiz/quiz_lecture5.md)      | [Quiz](quiz/quiz_lecture6.md)      
None      | None      | None      | None      | None      | None      
[YouTube](https://youtu.be/m_dle8vr3dU)      | [YouTube](https://youtu.be/B0wgOq_ECpI)      | [YouTube](https://youtu.be/AS2ubgAKW80)      | [YouTube](https://youtu.be/zDAMAdWk0DM)      | [YouTube](https://youtu.be/Jd86s2-8Zi0)      | [YouTube](https://youtu.be/XTalV2kUDPA)      


[Lecture 7](lecture7.md) | [Lecture 8](lecture8.md) | [Lecture 9](lecture9.md) | [Lecture 10](lecture10.md) | [Lecture 11](lecture11.md) | [Lecture 12](lecture12.md)
------------ | ------------- | ------------ | ------------- | ------------ | ------------- 
[Resources](resources/resource_lecture7.md) | [Resources](resources/resource_lecture8.md) | [Resources](resources/resource_lecture9.md) | [Resources](resources/resource_lecture10.md)  | [Resources](resources/resource_lecture11.md)  | None
[Quiz](quiz/quiz_lecture7.md)      | None      | None      | [Quiz](quiz/quiz_lecture10.md)       | None       | [Quiz](quiz/quiz_lecture12.md)
[Extra](extra/extra_lecture7.md)     | None      | None      | None       | None       | None
[YouTube](https://youtu.be/HoD9UJvOkPQ)     | [YouTube](https://youtu.be/b35LvlSZOuc)      | [YouTube](https://youtu.be/7FRiqSK6mBU)      | [YouTube](https://youtu.be/87KAvsfV73w)       | [YouTube](https://youtu.be/mrAoE4lFvas)       | [YouTube](https://youtu.be/IoTquSbTgoQ)

# Code Reference
The name of the codes/programs have been modified to recognize the pupose of the codes easily. The names used in the lectures can be found from the following table:

Lecture | GitHub Code Name | Original Name | GitHub Code Name | Original Name 
------- | --------- | ------------- | --------- | -------------
Lecture 1 | [fftw_1d.f95](fftw_1d.f95) | test.f95 | [fluid_1d.f95](fluid_1d.f95) | test_1.f95
Lecture 2 | [fourier_transform.f95](fourier_transform.f95) | fourier.f95 | [fluid_1d.f95](fluid_1d.f95) | test_1.f95
Lecture 3 | [fftw_2d.f95](fftw_2d.f95) | 2d_fft.f95
Lecture 4 | [poisson_2d.f95](poisson_2d.f95) | 2d_fft.f95
Lecture 5 | [fluid_2d.f95](fluid_2d.f95) | 2d_fft.f95
Lecture 6 | [fluid_2d.f95](fluid_2d.f95) | 2d_fft.f95
Lecture 7 |
Lecture 8 | [omp_do.f95](omp_do.f95) | omp.f95 | [omp_nested_do.f95](omp_nested_do.f95) | omp.f95
Lecture 9 | [mpi_do.f95](mpi_do.f95) | mpi.f95
Lecture 10| [mpi_nested_do.f95](mpi_nested_do.f95) | mpi.f95 | [hybrid_do.f95](hybrid_do.f95) | hybrid.f95
Lecture 11| [fftw_omp.f95](fftw_omp.f95) | 2d_fft.f95 | [fftw_mpi.f95](fftw_mpi.f95) | 2d_fft.f95
Lecture 12| [fluid_2d_with_tracers.f95](fluid_2d_with_tracers.f95) | 2d_fft.f95

# Code Description
## fourier_transform.f95: 
This code is a one dimensional serial Fourier transform solver. No external library is used in the code. It is for teaching / demonstration purpose only explicitly documenting the formulae of Fourier series from standard texts.
## fftw_1d.f95: 
This code takes a set of random numbers as input, takes their Fourier transform using one dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
## fftw_2d.f95: 
This code takes a set of random numbers as input of two dimensional array, takes their Fourier transform using two dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
## fftw_omp.f95: 
This is a program to teach how to use FFTW routine with OpenMP parallelization. This is a multi-core extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array in multiple cores, takes their Fourier transform using OpenMP parallel two dimensional FFTW library and then calculates the multi-core inverse Fourier transform, thus finally comparing the result with the input random numbers.
## fftw_mpi.f95: 
This is a program to teach how to use FFTW routine with MPI parallelization. This is a multi-node extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array, takes their Fourier transform using MPI parallel two dimensional FFTW library and then calculates the multi-node inverse Fourier transform, thus finally comparing the result with the input random numbers.
## poisson_2d.f95: 
This code evaluates the Poisson equation in two dimension using psedo-spectral method.
## fluid_1d.f95: 
This code solves one dimensional Burgers equation using pseudo-spectral method for spatial discretization and Adams-Bashforth algorithm for temporal update, with sin wave as initial condition. A shock appears at time t = 1.
## fluid_2d.f95: 
This code evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. The initial condition is taken as counter-streaming flows thereby akin to Kelvin-Helmholtz type instability.
## fluid_2d_with_tracers.f95: 
This code is an extension of the code fluid_2d.f95. Identical to the fluid_2d.f95, it also evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. However, this code has several passive tracer particles sprinkled in the fluid and the evolution of the passive tracer particles are followed using Cloud-In-Cell algorithm.
## omp_do.f95: 
This is a program to teach OpenMP do loops.
## omp_nested_do.f95: 
This is a program to teach OpenMP nested-do loops.
## mpi_do.f95: 
This is a program to teach MPI do loops.
## mpi_nested_do.f95: 
This is a program to teach MPI nested-do loops.
## hybrid_do.f95: 
This is a program to teach Hybrid (= OpenMP + MPI) nested-do loops.

# Commands to run the programs/codes
#### Serial Do Loop compilation
```console
gfortran <program_name.f95>; ./a.out
```
#### OpenMP Do Loop compilation
```console
gfortran -fopenmp <program_name.f95>; ./a.out
```
#### MPI Do Loop compilation
```console
mpif90 <program_name.f95>; mpirun -np 4 ./a.out
```
#### Hybrid Do Loop compilation
```console
mpif90 -fopenmp <program_name.f95>; mpirun -np 4 ./a.out
```
#### Serial FFTW compilation
```console
gfortran -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3 -lm; ./a.out
```
#### OpenMP FFTW compilation
```console
gfortran -fopenmp -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_omp -lfftw3 -lm; ./a.out
```
#### MPI FFTW compilation
```console
mpif90 -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_mpi -lfftw3 -lm; mpirun -np 4 ./a.out
```
#### Hybrid FFTW compilation
```console
mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun -quiet -np 4 ./a.out
```
#### Hybrid FFTW compilation and running with hyperthreads in macOS
```console
mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib <program_name.f95> -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun --use-hwthread-cpus ./a.out
```
