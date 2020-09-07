# Fluid teaching
a repository for teaching algorithms for fluid codes!

Lecture 1: 
Lecture 2: 
Lecture 3: 
Lecture 4: 
Lecture 5: 
Lecture 6: 
Lecture 7: Revision
Lecture 8: 
Lecture 9: 
Lecture 10: 
Lecture 11: 
Lecture 12: 

=== Quiz ===

Lecture 1: No Quiz

Lecture 2: No Quiz

Lecture 3: 
Quiz - 1: Can you modify the 1D Burgers code, and reproduce electron-plasma-oscillation? [Hint: Look at the plasma oscillation section of Davidson's book]
Quiz - 2: And eventually increase the amplitude of perturbation and see how does the plasma frequency changes as nonlinearity enters via large amplitude perturbation. Remember, our code successfully passed one of the most difficult tests of nonlinearity - the shock problem!!!

Lecture 4: No Quiz
Well, I just identified the error.
Look at Line: 347 in the Subroutine AdamsBashforth.
    omegak_new(i,j) = omegak(Nh,Ny) + ..........
should be replaced as
    omegak_new(i,j) = omegak(i,j) + ...........
You can now replace it in the shared version and try to check!

Lecture 5: https://en.wikipedia.org/wiki/Ackermann_function

:<math> 
\begin{array}{lcl}
\operatorname{A}(0, n) & = & n + 1 \\
\operatorname{A}(m+1, 0) & = & \operatorname{A}(m, 1) \\
\operatorname{A}(m+1, n+1) & = & \operatorname{A}(m, A(m+1, n))
\end{array}
</math>

Evaluate A(3,11) and send your answer to me (each of you individually)! 
Do NOT use 'function' call.

Lecture 6: 
Quiz - 4: Since all of you now know, how to take 1D and 2D Fourier transforms, can you now create a 3D array and take its Fourier transform and then inverse Fourier transform and check whether you get back the input 3D array?

And, regarding Aliasing and de-aliasing methods, I have explained in short in the next lecture. But for a more detailed overview you can look at the following link:
https://en.wikipedia.org/wiki/Aliasing
https://www.astro.auth.gr/~vlahos/GravitoplasmaWS1/pseudo-spectral_2.pdf

Resources:
https://arxiv.org/abs/1711.10865
https://doi.org/10.1073/pnas.1509304112
https://doi.org/10.1103/PhysRevLett.75.2486
https://doi.org/10.1017/S0022112061000378
https://doi.org/10.1088/1742-6596/1548/1/012037

Lecture 7: Revision
Attaching few papers on simulating bounded flows via pseudo-spectral method. May be we can have a bit of discussion sometime later on such simulations, since flows in bounded domains are the most natural systems. Also, you guys can try to implement the technique in the 2D code that all of us together wrote few weeks ago! You can also search independently about such methods. One helpful google-search key may be - "Volume penalization method".:slightly_smiling_face:
One other option to simulate boundaries / bounded flows, within pseudo-spectral scheme is - instead of using the Fourier transform, one can use sine-Fourier or cosine-Fourier transforms while taking the derivatives. Such libraries are easily available within the FFTW architecture but I had a hard time in the implementation process. May be, young blood can rejuvenate the endeavor!!!
http://dx.doi.org/10.1016/j.jcp.2014.05.038
https://doi.org/10.1016/j.cpc.2010.05.019
https://doi.org/10.1016/j.compfluid.2004.09.006

You can get some funny videos on turbulence (some of which you can now simulate now using the 2dfft.f95 code) here.
https://www.youtube.com/watch?v=5zI9sG3pjVU
In the middle of the week, I thought of sharing some 'funny' movies!
https://w3.pppl.gov/~hammett/viz/viz.html
You can watch the first two movies and find out what is/are the primary instability(s) occurring within the DIII-D tokamak!
(DIII-D is just a name of a tokamak, somewhere in the west coast of USA, in case you have not heard about it earlier)
And this takes me to Quiz - 6:
Can you write down, what extra things we need to add in our 2dfft.f95 code, to simulate such a 'real-life' plasma?
So thought of sharing the link with you all.
https://www.youtube.com/playlist?list=PLYwpaL_SFmcA1eJbqwvjKgsnT321hXRGx
Lecture 8: 
Lecture 9: 
Lecture 10: 
https://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-loop.html
http://users.metu.edu.tr/csert/me582/ME582%20Ch%2001.pdf
Dear all,
I added some comments at the end of the program that we just wrote today.
If you can run the program, it will print some comments those comments will give you some hint about, why for 2 nodes, it did not give us correct result.
Also, further at the end, it will print your Quiz - 8!
And finally what if, I send you a file and ask, what does this program do?
To compile and run this file, you may use the following command:
mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib 2d_fft.f95 -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun -quiet -np 3 ./a.out
Quiz - 9: Can you turn the code [3d_hybrid_poisson.f95] into a 3d hybrid Poisson solver?
If you can solve Quiz - 9, this course is over!!!
Lecture 11: 
Here is a nice link from where, I learnt OpenMP first. May be you guys also will like:
https://chryswoods.com/beginning_openmp/index.html
And here is his general course link:
https://chryswoods.com/main/courses.html
Lecture 12: 
Quiz - ??: The 2D algorithm I described today for the tracer particles (interpolating from grid to particle position) had 4 if-conditions . Now check out the 'big' expression (https://arxiv.org/pdf/1810.12707.pdf, Page No 7, Section VI, A, at the bottom part of the page) and tell me how many if-conditions  you have to write for this case?


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
