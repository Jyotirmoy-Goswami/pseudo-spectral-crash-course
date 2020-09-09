Fluid teaching
===============================================
A repository of teaching numerical algorithms for solving fluid equations. 

We have described pseudo-spectral method for spatial discretization. The equations solved in this lecture series are two dimensional incompressible fluid equations in vorticity (<img src="https://render.githubusercontent.com/render/math?math=\omega">) - stream function (<img src="https://render.githubusercontent.com/render/math?math=\psi">) formalism.


<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial \vec{\omega}}{\partial t} %2B \vec{u} \cdot \vec{\nabla} \vec{\omega} = \nu \nabla^2 \vec{\omega}">
<img src="https://render.githubusercontent.com/render/math?math=\vec{\omega} = \vec{\nabla} \times \vec{u}">
<img src="https://render.githubusercontent.com/render/math?math=\nabla^2 \psi = - \omega ">


The wave-form is assumed to be <img src="https://render.githubusercontent.com/render/math?math=\exp(i(kx-wt))">. Hence derivative w.r.t. x gives ik only (and NOT -ik).

A series of lectures are available on youtube based on the materials archived in this repository. Click on the link provided below to access the lectures.


[Webinar on "Pseudo Spectral Method"](https://www.youtube.com/playlist?list=PLbX_ZyxeXxSJXnIAnkhhAsIAV-Ld0Awsu)

[![Webinar on "Pseudo Spectral Method"](http://img.youtube.com/vi/m_dle8vr3dU/0.jpg)](https://www.youtube.com/embed/videoseries?list=PLbX_ZyxeXxSJXnIAnkhhAsIAV-Ld0Awsu)
# Lecture description
## Lecture 1: 
To make a post-programming fun of my coding, you guys can go to Time = 1:41:33 of the first video link and see how the error entered into our code due to copy-paste from the previous line, which kept us bothering!:joy:

I was laughing, watching myself making the mistake. It is really funny and amusing. This is my life's first experience when, I can 'post-see' myself making a mistake.:laughing:
## Lecture 2: 
Yet to be updated
## Lecture 3: 
Yet to be updated
## Lecture 4: 
Yet to be updated
## Lecture 5: 
Yet to be updated
## Lecture 6: 
Yet to be updated
## Lecture 7 (Revision): 
Yet to be updated
## Lecture 8: 
Yet to be updated
## Lecture 9: 
Yet to be updated
## Lecture 10: 
Yet to be updated
## Lecture 11: 
Yet to be updated
## Lecture 12: 
Yet to be updated


# Resources and tips related to each lecture

## Lecture 1: 
No resources Available
## Lecture 2: 
As mentioned, I add the link of Blackboard lectures by JKB and Rama Govindarajan.

[Lecture 1 by JKB](https://www.youtube.com/watch?v=0JMeOwgQT-k&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=4)

[Lecture 1 by Rama Govindarajan](https://www.youtube.com/watch?v=-zbwHOXiLzc&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=6)

Follow-up lectures (Lecture 2,3,4 etc.) should be found in the upper-right-hand corner.

And, and, and...

Here you can find, how I am cheating you in every lecture.:wink: [And hopefully, in the next lecture also](https://www.youtube.com/watch?v=EOYJc2Vuju0&index=54&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&hd=1):shushing_face:

Lastly, kudos to all of you, for your immense patience in my lousy discussions. Each of you really deserve at-least one candy!:yum:

So, [here is yours!](https://www.youtube.com/watch?v=tFtpM-Evo90&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=37)

Happy coding!
## Lecture 3: 
No resources Available
## Lecture 4: 
There was an error in the code we were writing today :smiley:. Well, I just identified the error. 

Look at **Line: 347** in the Subroutine *AdamsBashforth*.
```
    omegak_new(i,j) = omegak(Nh,Ny) + ..........
```
should be replaced as
```
    omegak_new(i,j) = omegak(i,j) + ...........
```
## Lecture 5: 
No resources Available
## Lecture 6: 

Regarding Aliasing and de-aliasing methods, I will explain in short in the next lecture (**Lecture 7**). But for a more detailed overview you can look at the following link:
* https://en.wikipedia.org/wiki/Aliasing
* https://www.astro.auth.gr/~vlahos/GravitoplasmaWS1/pseudo-spectral_2.pdf

### Resources:
* https://arxiv.org/abs/1711.10865
* https://doi.org/10.1073/pnas.1509304112
* https://doi.org/10.1103/PhysRevLett.75.2486
* https://doi.org/10.1017/S0022112061000378
* https://doi.org/10.1088/1742-6596/1548/1/012037

## Lecture 7: 
Following are few papers on simulating bounded flows via pseudo-spectral method. May be we can have a bit of discussion sometime later on such simulations, since flows in bounded domains are the most natural systems. Also, you guys can try to implement the technique in the 2D code that all of us together wrote few weeks ago! You can also search independently about such methods. One helpful google-search key may be - "Volume penalization method".

One other option to simulate boundaries / bounded flows, within pseudo-spectral scheme is - instead of using the Fourier transform, one can use sine-Fourier or cosine-Fourier transforms while taking the derivatives. Such libraries are easily available within the FFTW architecture but I had a hard time in the implementation process. May be, young blood can rejuvenate the endeavor!!!
### Resources:
* http://dx.doi.org/10.1016/j.jcp.2014.05.038
* https://doi.org/10.1016/j.cpc.2010.05.019
* https://doi.org/10.1016/j.compfluid.2004.09.006

You can get some funny videos on turbulence (some of which you can now simulate now using the 2dfft.f95 code) here.
https://www.youtube.com/watch?v=5zI9sG3pjVU

And here is one more link (http://www.lcs-fast.com/fifth_order/) that shows some spectral simulation of air flows around a racing-car. But note that these are not pseudo-spectral. The basic idea is same though. The only difference is, instead of Fourier basis, they have used Legendre or Chebyshev basis (and sometimes Hermite also!) for better accuracy! If you are feeling confused, revise our first discussion and specifically look at the first slide!
##### Extras
In the middle of the week, I thought of sharing some 'funny' movies!
https://w3.pppl.gov/~hammett/viz/viz.html
You can watch the first two movies and find out what is/are the primary instability(s) occurring within the DIII-D tokamak!
(DIII-D is just a name of a tokamak, somewhere in the west coast of USA, in case you have not heard about it earlier)

And this takes me to [Quiz - 6](#quiz---6):

## Lecture 8: 
Need inspiration? watch this,
[High Performance Computing /Parallel Computing](https://www.youtube.com/playlist?list=PLYwpaL_SFmcA1eJbqwvjKgsnT321hXRGx)
## Lecture 9: 

Here is a question cum hint for you to proceed for [Quiz - 1](#quiz---1):
What is the dimension of x in the expression exp(x)?

And here is a link of Davidson's book [Methods in Nonlinear Plasma Theory](https://books.google.no/books?id=8iW0MDOVr0oC&lpg=PP1&hl=no&pg=PP1#v=onepage&q&f=false). You may go to Chapter 3, Section 3.1, Page 33-34 and read-up a little.

## Lecture 10: 


* https://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-loop.html
* http://users.metu.edu.tr/csert/me582/ME582%20Ch%2001.pdf

## Lecture 11: 

Here is a nice link from where, I learnt OpenMP first. May be you guys also will like:
* https://chryswoods.com/beginning_openmp/index.html

And here is his general course link:
* https://chryswoods.com/main/courses.html

## Lecture 12: 


# Quizzes 

## Lecture 1: 
No Quiz
## Lecture 2: 
No Quiz
## Lecture 3: 
### Quiz - 1: 
Can you modify the 1D Burgers code, and reproduce electron-plasma-oscillation? (Hint: Look at the plasma oscillation section of Davidson's book)
### Quiz - 2: 
And eventually increase the amplitude of perturbation and see how does the plasma frequency changes as nonlinearity enters via large amplitude perturbation. Remember, our code successfully passed one of the most difficult tests of nonlinearity - the shock problem!!!
## Lecture 4: 
No Quiz
## Lecture 5: 
### Quiz - 3:
[Ackermann function](https://en.wikipedia.org/wiki/Ackermann_function)

<img src="https://render.githubusercontent.com/render/math?math=A(0,n) = (n%2B1)">

<img src="https://render.githubusercontent.com/render/math?math=A(m%2B1,0) = A(m,1)">

<img src="https://render.githubusercontent.com/render/math?math=A(m%2B1,n%2B1) = A(m,A(m%2B1,n))">

Evaluate A(3,11) and do NOT use *function* call.

## Lecture 6: 
### Quiz - 4: 
Since all of you now know, how to take 1D and 2D Fourier transforms, can you now create a 3D array and take its Fourier transform and then inverse Fourier transform and check whether you get back the input 3D array?

## Lecture 7:
No Quiz
### Quiz - 5: 
Solution: derivative_3d_recipe.md
### Quiz - 6: 
Can you write down, what extra things we need to add in our 2dfft.f95 code, to simulate such a 'real-life' plasma in these funny movies described in [Extras](#extras)?

## Lecture 8: 
No Quiz
## Lecture 9: 
No Quiz

## Lecture 10: 

I added some comments at the end of the program that we just wrote today. If you can run the program, it will print some comments those comments will give you some hint about, why for 2 nodes, it did not give us correct result.
Also, further at the end, it will print your Quiz - 8!

And finally what if, I send you a file and ask, what does [this program](https://github.com/RupakMukherjee/fluid_teaching/blob/master/3d_hybrid_poisson.f95) do?
To compile and run this file, you may use the following command:
```console
mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib 3d_hybrid_poisson.f95 -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun -quiet -np 3 ./a.out
```
### Quiz - 9: 
Can you turn the code [3d_hybrid_poisson.f95] into a 3d hybrid Poisson solver?
If you can solve Quiz - 9, this course is over!!!

## Lecture 11: 
No Quiz

## Lecture 12: 
### Quiz - 10: 
The 2D algorithm I described today for the tracer particles (interpolating from grid to particle position) had 4 if-conditions . Now check out the 'big' expression (https://arxiv.org/pdf/1810.12707.pdf, Page No 7, Section VI, A, at the bottom part of the page) and tell me how many if-conditions  you have to write for this case?

# Code Reference
The name of the codes/programs have been modified to recognize the pupose of the codes easily. The names used in the lectures can be found from the following table:

Lecture | Code Name | Original Name | Code Name | Original Name 
------- | --------- | ------------- | --------- | -------------
Lecture 1 | fftw_1d.f95 | test.f95 | burgulence.f95 | test_1.f95
Lecture 2 | fourier_transform.f95 | fourier.f95 | burgulence.f95 | test_1.f95
Lecture 3 | fftw_2d.f95 | 2d_fft.f95
Lecture 4 | poisson_2d.f95 | 2d_fft.f95
Lecture 5 | fluid_2d.f95 | 2d_fft.f95
Lecture 6 | fluid_2d.f95 | 2d_fft.f95
Lecture 7 |
Lecture 8 | omp_do.f95 | omp.f95 | omp_nested_do.f95 | omp.f95
Lecture 9 | mpi_do.f95 | mpi.f95
Lecture 10| mpi_nested_do.f95 | mpi.f95 | hybrid_do.f95 | hybrid.f95
Lecture 11| fftw_omp.f95 | 2d_fft.f95 | fftw_mpi.f95 | 2d_fft.f95
Lecture 12| fluid_2d_with_tracers.f95 | 2d_fft.f95

# Code Description
## fftw_1d.f95: 
This code takes a set of random numbers as input, takes their Fourier transform using one dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
## burgulence.f95: 
This code solves one dimensional Burgers equation using pseudo-spectral method for spatial discretization and Adams-Bashforth algorithm for temporal update, with sin wave as initial condition. A shock appears at time t = 1.
## fourier_transform.f95: 
This code is a one dimensional serial Fourier transform solver. No external library is used in the code. It is for teaching / demonstration purpose only explicitly documenting the formulae of Fourier series from standard texts.
## fftw_2d.f95: 
This code takes a set of random numbers as input of two dimensional array, takes their Fourier transform using two dimensional FFTW library and then calculates the inverse Fourier transform, thus finally comparing the result with the input random numbers.
## poisson_2d.f95: 
This code evaluates the Poisson equation in two dimension using psedo-spectral method.
## fluid_2d.f95: 
This code evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. The initial condition is taken as counter-streaming flows thereby akin to Kelvin-Helmholtz type instability.
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
## fftw_omp.f95: 
This is a program to teach how to use FFTW routine with OpenMP parallelization. This is a multi-core extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array in multiple cores, takes their Fourier transform using OpenMP parallel two dimensional FFTW library and then calculates the multi-core inverse Fourier transform, thus finally comparing the result with the input random numbers.
## fftw_mpi.f95: 
This is a program to teach how to use FFTW routine with MPI parallelization. This is a multi-node extension of fftw_2d.f95. Hence this code also takes a set of random numbers as input of two dimensional array, takes their Fourier transform using MPI parallel two dimensional FFTW library and then calculates the multi-node inverse Fourier transform, thus finally comparing the result with the input random numbers.
## fluid_2d_with_tracers.f95: 
This code is an extension of the code fluid_2d.f95. Identical to the fluid_2d.f95, it also evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. However, this code has several passive tracer particles sprinkled in the fluid and the evolution of the passive tracer particles are followed using Cloud-In-Cell algorithm.

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
