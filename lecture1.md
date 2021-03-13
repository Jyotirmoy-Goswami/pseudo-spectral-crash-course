# Lecture 1: 
In the first part of the lecture, I described some basic idea about different spectral methods, concept of basis function, projection of a function into different basis, accuracy of spectral methods over finite difference schemes, idea of spectral convergence. Further I chose one of the many available complete basis function (i.e. Fourier basis) and talked about pseudo-spectral method, specifically showing one example, named "Burgers equation". The analytical solution of this equation is [well documented](https://arxiv.org/pdf/nlin/0012033.pdf).

In the last part of the lecture, I tried to numerically solve the one dimensional Burgers equation using pseudo-spectral method (for Fourier transforms I used [FFTW library](http://www.fftw.org/)) which failed.:joy: Before going into the pseudo-spectral method, I also described how to use the one dimensional Fourier transform using FFTW library.

To make a post-programming fun of my coding, you can go to [Time = 1:41:33](https://youtu.be/m_dle8vr3dU?t=6093) of the first video link and see how the error entered into our code due to copy-paste from the previous line, which kept us bothering!:joy:

I was laughing, watching myself making the mistake. It is really funny and amusing. This is my life's first experience when, I can 'post-see' myself making a mistake.:laughing:

Anyway, to take a look at what went wrong into our code, you can always check-out the ([Part-2](https://www.youtube.com/watch?v=csuAlbc-4cg&list=PLbX_ZyxeXxSJXnIAnkhhAsIAV-Ld0Awsu)) of the Lecture 1!

