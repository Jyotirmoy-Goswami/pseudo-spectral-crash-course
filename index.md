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

[Lecture 1](lecture1.md) | [Lecture 2](lecture2.md) | [Lecture 3](lecture3.md) | [Lecture 4](lecture4.md) | [Lecture 5](lecture5.md) | [Lecture 6](lecture6.md) | [Lecture 7](lecture7.md) | [Lecture 8](lecture8.md) | [Lecture 9](lecture9.md) | [Lecture 10](lecture10.md) | [Lecture 11](lecture11.md) | [Lecture 12](lecture12.md)
------------ | ------------- | ------------ | ------------- | ------------ | ------------- | ------------ | ------------- | ------------ | ------------- | ------------ | ------------- 
[Resources](resources/resource_lecture1.md) | [Resources] | None      | None      | None      | [Resources] | [Resources] | [Resources] | [Resources] | [Resources]  | [Resources]  | None
None      | None      | [Quiz]      | None      | [Quiz]      | [Quiz]      | [Quiz]      | None      | None      | [Quiz]       | None       | [Quiz]
None      | None      | None      | None      | None      | None      | [Extra]     | None      | None      | None       | None       | None

          
## Lecture 1: 
In the first part of the lecture, I described some basic idea about different spectral methods, concept of basis function, projection of a function into different basis, accuracy of spectral methods over finite difference schemes, idea of spectral convergence. Further I chose one of the many available complete basis function (i.e. Fourier basis) and talked about pseudo-spectral method, specifically showing one example, named "Burgers equation". The analytical solution of this equation is [well documented](https://arxiv.org/pdf/nlin/0012033.pdf).

In the last part of the lecture, I tried to numerically solve the one dimensional Burgers equation using pseudo-spectral method (for Fourier transforms I used [FFTW library](http://www.fftw.org/)) which failed.:joy: Before going into the pseudo-spectral method, I also described how to use the one dimensional Fourier transform using FFTW library.

To make a post-programming fun of my coding, you can go to [Time = 1:41:33](https://youtu.be/m_dle8vr3dU?t=6093) of the first video link and see how the error entered into our code due to copy-paste from the previous line, which kept us bothering!:joy:

I was laughing, watching myself making the mistake. It is really funny and amusing. This is my life's first experience when, I can 'post-see' myself making a mistake.:laughing:

Anyway, to take a look at what went wrong into our code, you can always check-out the ([Part-2](https://youtu.be/VDcAC_IW8s4)) of the Lecture 1!



<details>
<summary>Resources</summary>  
Instruction for installing FFTW library:

Download the latest version from [here](http://www.fftw.org/download.html)

```console
$ cd path_to_file
$ sudo ./configure --enable-threads --enable-openmp --with-g77-wrappers
$ sudo make
$ sudo make install
```
</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/m_dle8vr3dU)


## Lecture 2: 
In this lecture, I first delineated some basic properties of the Burgers equation and compared them with our freshly brewed "[fluid_1d.f95](fluid_1d.f95)" code (the one that we developed in the first lecture) and then described some basic idea about discontinuity (which will appear in the form of 'shock' in later part of this course) and how to capture that using our simple numerical techniques. Next I compare the algorithm mentioned in my slide with the code we developed in the last lecture. And then got drifted into turbulence, scaling, energy cascades etc. which I did not intend to cover.

<details>
<summary>Resources</summary>  
    
As mentioned, I add the link of Blackboard lectures by Jayanta K. Bhattacharjee(JKB) and Rama Govindarajan.

[Lecture 1 by JKB](https://www.youtube.com/watch?v=0JMeOwgQT-k&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=4)

[Lecture 1 by Rama Govindarajan](https://www.youtube.com/watch?v=-zbwHOXiLzc&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=6)

Follow-up lectures (Lecture 2,3,4 etc.) should be found in the upper-right-hand corner.

And, and, and...

Here you can find, how I am cheating you in every lecture.:wink: [And hopefully, in the next lecture also](https://www.youtube.com/watch?v=EOYJc2Vuju0&index=54&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&hd=1):shushing_face:

Lastly, kudos to all of you, for your immense patience in my lousy discussions. Each of you really deserve at-least one candy!:yum:

So, [here is yours!](https://www.youtube.com/watch?v=tFtpM-Evo90&list=PL04QVxpjcnjhcA2iryoRvU86o939OrZSo&index=37)

Happy coding!

</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/B0wgOq_ECpI)

## Lecture 3: 
Here I talked about the two dimensional form of the Navier-Stokes equation, its component-by-component decomposition into Cartesian basis, the vorticity - stream-function formalism for incompressible fluid and finally some tips on how to use the multi-dimensional FFTW library. As a prelude to pseudo-spectrally solve the two dimensional incompressible Navier-Stokes equation using multi-dimensional FFTW library, I first took a set of uniform random numbers as an input of a 2D array, took its Fourier transform and then took inverse Fourier transform to come back to the real space. If I get back the input array, it shows that we have learned using the multi-dimensional FFTW library perfectly! In the last part, I described a bit about the two dimensional Poisson solver using this freshly brewed code! The numerical example is solved in the next lecture.

<details>
<summary>Quizzes</summary> 
    
### Quiz - 1: 
Can you modify the [1D Burgers code](burgulence.f95), and reproduce electron-plasma-oscillation? 

(Hint: Look at the plasma oscillation section of Davidson's book)

### Quiz - 2: 
And eventually increase the amplitude of perturbation and see how does the plasma frequency changes as nonlinearity enters via large amplitude perturbation. Remember, our code successfully passed one of the most difficult tests of nonlinearity - the shock problem!!!

</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/AS2ubgAKW80)

## Lecture 4: 
As promised, I solved the two dimensional Poisson solver using multi-dimensional FFTW library, but before that, I revised the numerical subtleties of multi-dimensional FFTW library that I described in the last lecture. The revision took almost half of the lecture and in the rest of the part I showed the numerical implementation as an extension of the program written in the last lecture. At the fag end I talked about some seminal benchmarking papers for two dimensional fluid codes!

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/zDAMAdWk0DM)

## Lecture 5: 
This is probably the most crucial lecture. Starting from last [lecture's code](poisson_2d.f95), I developed a two dimensional incompressible fluid solver. In the first part I described the basic equaitons to be solved in the vorticity - stream-function formalism at length and then developed the code. Unfortunately, the code had a bug which I later identified and pointed out in the beginning of the next lecture.

<details>
<summary>Quizzes</summary> 
    
### Quiz - 3:
[Ackermann function](https://en.wikipedia.org/wiki/Ackermann_function)


* <img src="https://render.githubusercontent.com/render/math?math=A(0,n) = (n%2B1)">


* <img src="https://render.githubusercontent.com/render/math?math=A(m%2B1,0) = A(m,1)">


* <img src="https://render.githubusercontent.com/render/math?math=A(m%2B1,n%2B1) = A(m,A(m%2B1,n))">


Evaluate <img src="https://render.githubusercontent.com/render/math?math=A(3,11)"> and do NOT use *function* call.

</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/Jd86s2-8Zi0)

## Lecture 6: 
There was an error in the code we were writing in the last lecture :smiley:. Well, I identified the error. 

Look at **Line: 347** in the Subroutine *AdamsBashforth*.
```
    omegak_new(i,j) = omegak(Nh,Ny) + ..........
```
should be replaced as
```
    omegak_new(i,j) = omegak(i,j) + ...........
```
Afterwards, I talked about the basic idea of Kelvin-Helmholtz instability, how to capture it numerically (there are some subtleties here). Finally our baby-code performed well to reproduce the analytical growth rate calculated by P G Drazin and thus passed one very critical test! Then I mentioned some papers, talked about few interesting open problems that can be 'easily' attacked using our simple baby code and finally sketched some aspects of parallel computing that we will be covering in the rest of the course.


<details>
<summary>Resources</summary> 
    
Regarding Aliasing and de-aliasing methods, I will explain in short in the next lecture (**Lecture 7**). But for a more detailed overview you can look at the following link:
* https://en.wikipedia.org/wiki/Aliasing
* https://www.astro.auth.gr/~vlahos/GravitoplasmaWS1/pseudo-spectral_2.pdf

### Resources:
* https://arxiv.org/abs/1711.10865
* https://doi.org/10.1073/pnas.1509304112
* https://doi.org/10.1103/PhysRevLett.75.2486
* https://doi.org/10.1017/S0022112061000378
* https://doi.org/10.1088/1742-6596/1548/1/012037

</details>


<details>
<summary>Quizzes</summary> 
    
### Quiz - 4: 
Since all of you now know, how to take 1D and 2D Fourier transforms, can you now create a 3D array and take its Fourier transform and then inverse Fourier transform and check whether you get back the input 3D array?

Hint:

```fortran
do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    do k = 1,Nz/2
      kz = 2.0d0*pi*dfloat(k-1)/Lz 
      i_kx_omegak(i,j,k) = (0.0d0,1.0d0) * kx * omegak(i,j,k)
      i_ky_omegak(i,j,k) = (0.0d0,1.0d0) * ky * omegak(i,j,k)
      i_kz_omegak(i,j,k) = (0.0d0,1.0d0) * kz * omegak(i,j,k)      
    enddo
    do k = Nz/2+1,Nz
      kz = 2.0d0*pi*(dfloat(k-1)-Nz)/Lz
      i_kx_omegak(i,j,k) = (0.0d0,1.0d0) * kx * omegak(i,j,k)
      i_ky_omegak(i,j,k) = (0.0d0,1.0d0) * ky * omegak(i,j,k)
      i_kz_omegak(i,j,k) = (0.0d0,1.0d0) * kz * omegak(i,j,k)             
    enddo    
  do j = Ny/2+1,Ny
     ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    do k = 1,Nz/2
      kz = 2.0d0*pi*dfloat(k-1)/Lz 
      i_kx_omegak(i,j,k) = (0.0d0,1.0d0) * kx * omegak(i,j,k)
      i_ky_omegak(i,j,k) = (0.0d0,1.0d0) * ky * omegak(i,j,k)
      i_kz_omegak(i,j,k) = (0.0d0,1.0d0) * kz * omegak(i,j,k)      
    enddo
    do k =  Nz/2+1,Nz
      kz = 2.0d0*pi*(dfloat(k-1)-Nz)/Lz
      i_kx_omegak(i,j,k) = (0.0d0,1.0d0) * kx * omegak(i,j,k)
      i_ky_omegak(i,j,k) = (0.0d0,1.0d0) * ky * omegak(i,j,k)
      i_kz_omegak(i,j,k) = (0.0d0,1.0d0) * kz * omegak(i,j,k)             
    enddo     
  enddo
enddo
```

</details>


Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/XTalV2kUDPA)


## Lecture 7 (Revision): 
I got a feedback that I am going too fast. Hence we decided to pause for a bit, put one full lecture for revision and then proceed for parallelization. So in this lecture, I talked again about the numerical aspects of capturing Kelvin-Helmholtz instability, analytical growth rate, numerical comparison, delineated a bit in to aliazing errors and went back into the lecture slide of the first lecture, compared our previous codes line-by-line with the algorithm described in the lecture slide and finally for the first time touched the aspects of three dimensional arrays.

##### Extras
In the middle of the week, I thought of sharing some 'funny' movies!
https://w3.pppl.gov/~hammett/viz/viz.html
You can watch the first two movies and find out what is/are the primary instability(s) occurring within the DIII-D tokamak!
(DIII-D is just a name of a tokamak, somewhere in the west coast of USA, in case you have not heard about it earlier)

And this takes me to Quizzes of this lecture.

<details>
<summary>Resources</summary>  
    
Following are few papers on simulating bounded flows via pseudo-spectral method. May be we can have a bit of discussion sometime later on such simulations, since flows in bounded domains are the most natural systems. Also, you guys can try to implement the technique in the 2D code that all of us together wrote few weeks ago! You can also search independently about such methods. One helpful google-search key may be - "Volume penalization method".

One other option to simulate boundaries / bounded flows, within pseudo-spectral scheme is - instead of using the Fourier transform, one can use sine-Fourier or cosine-Fourier transforms while taking the derivatives. Such libraries are easily available within the FFTW architecture but I had a hard time in the implementation process. May be, young blood can rejuvenate the endeavor!!!
### Resources:
* http://dx.doi.org/10.1016/j.jcp.2014.05.038
* https://doi.org/10.1016/j.cpc.2010.05.019
* https://doi.org/10.1016/j.compfluid.2004.09.006

You can get some funny videos on turbulence (some of which you can now simulate now using the 2dfft.f95 code) here.
https://www.youtube.com/watch?v=5zI9sG3pjVU

And here is one more link (http://www.lcs-fast.com/fifth_order/) that shows some spectral simulation of air flows around a racing-car. But note that these are not pseudo-spectral. The basic idea is same though. The only difference is, instead of Fourier basis, they have used Legendre or Chebyshev basis (and sometimes Hermite also!) for better accuracy! If you are feeling confused, revise our first discussion and specifically look at the first slide!

</details>


<details>
<summary>Quizzes</summary> 
    
### Quiz - 5: 
Can you write down, what extra things we need to add in our 2dfft.f95 code, to simulate such a 'real-life' plasma in these funny movies described in [Extras](#extras)?

</details>


Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/HoD9UJvOkPQ)

## Lecture 8: 
Here I started with some 'funny' movies (both hydrodynamics and plasmas) and don't know how, got drifted towards my favourite stereographic projection and topology while talking about periodic boundary condition. Finally I could hold myself and came back to simulation discussions on bounded flows using spectral methods. In the second half of the lecture, I talked about basic ideas about multi-core and multi-node parallelization and its simple numerical implementations, explicitly showing how to parallelize "do-loops" for OpenMP and MPI architechtures.

<details>
    
<summary>Resources</summary>  
Need inspiration? watch this,

[High Performance Computing /Parallel Computing](https://www.youtube.com/playlist?list=PLYwpaL_SFmcA1eJbqwvjKgsnT321hXRGx)

</details>


Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/b35LvlSZOuc)


## Lecture 9: 
In this lecture I talked about MPI architechture in very detail explaining the bits as elaborately as possible. At the end, I just had re-done the single "do-loop" example that I outlined in the last lecture. 

<details>
    
<summary>Resources</summary> 

Here is a question cum hint for you to proceed for [Quiz - 1](#quiz---1):
What is the dimension of x in the expression exp(x)?

And here is a link of Davidson's book [Methods in Nonlinear Plasma Theory](https://books.google.no/books?id=8iW0MDOVr0oC&lpg=PP1&hl=no&pg=PP1#v=onepage&q&f=false). You may go to Chapter 3, Section 3.1, Page 33-34 and read-up a little.

</details>


Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/7FRiqSK6mBU)

## Lecture 10: 
Here I went to 2D arrays! Basically I extended the "do-loop" example I described in the last lecture for two dimensional arrays. This extension is crucial since we will need it when you will be parallelizing the code we wrote up in Lecture 5 and 6.

<details>
    
<summary>Resources</summary> 

Few resources to read:
* https://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-loop.html
* http://users.metu.edu.tr/csert/me582/ME582%20Ch%2001.pdf

</details>


<details>
<summary>Quizzes</summary> 
I added some comments at the end of the program that we just wrote today. If you can run the program, it will print some comments those comments will give you some hint about, why for 2 nodes, it did not give us correct result.
    
Also, further at the end, it will print your Quiz - 8!

And finally what if, I send you a file and ask, what does [this program](https://github.com/RupakMukherjee/fluid_teaching/blob/master/poisson_3d_hybrid.f95) do?
To compile and run this file, you may use the following command:

```console
$ mpif90 -fopenmp -I/usr/local/include -L/usr/local/lib poisson_3d_hybrid.f95 -lfftw3_mpi -lfftw3_omp -lfftw3 -lm; mpirun -quiet -np 3 ./a.out
```

### Quiz - 6: 
Can you turn the code "[hybrid_do.f95]"(hybrid_do.f95) into a 3d hybrid Poisson solver?
If you can solve Quiz - 9, this course is over!!!

</details>


Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/87KAvsfV73w)


## Lecture 11: 
This is the most important lecture for you, if you want to get our code parallelized that we developed in Lecture 5 and 6. When I say 'parallelize', I mean both multi-core (OpenMP) and multi-node (MPI). After talking about how to parallelize nested "do-loops" in the last lecture, here I showed how to use the parallelized FFTW library. Thus at the end of the lecture, all of you become potentially capable to write an OpenMP and MPI parallel two dimensional incompressible Navier-Stokes equation solver! Start with the code we developed in [Lecture 5](#lecture-5) and [6](#lecture-6). Make all the "do-loops" parallel using the prescription in [Lecture 9](#lecture-9) and replace all the FFTW calls using the prescription in [Lecture 10](#lecture-10). You are all set!!!

<details>
    
<summary>Resources</summary> 

Here is a nice link from where, I learnt OpenMP first. May be you guys also will like:
* https://chryswoods.com/beginning_openmp/index.html

And here is his general course link:
* https://chryswoods.com/main/courses.html

</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/mrAoE4lFvas)


## Lecture 12: 
The first part of this lecture extends over a set of previous lecture series on molecular-dynamics simulation and then smoothly connects to "smooth particle hydrodynamics" and enters into the pseudo-spectral code we developed in [Lecture 5](#lecture-5) and [6](#lecture-6). Finally I talked about adding some passive tracer particles delineating some crude and very basic idea about particle-in-cell algorithm and particle pusher schemes we used in the molecular-dynamics simulation series lecture. In between the lecture I also have stressed the need of the implementation of this scheme and its potential applications!


<details>
    
<summary>Quizzes</summary> 

### Quiz - 7: 

The 2D algorithm I described today for the tracer particles (interpolating from grid to particle position) had 4 if-conditions . Now check out the 'big' expression (https://arxiv.org/pdf/1810.12707.pdf, Page No 7, Section VI, A, at the bottom part of the page) and tell me how many if-conditions  you have to write for this case?

</details>

Lecture Video: 

[<img src="yt_logo_rgb_light.png" width="100">](https://youtu.be/IoTquSbTgoQ)


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
