# Resource (Lecture 7)

Following are few papers on simulating bounded flows via pseudo-spectral method. May be we can have a bit of discussion sometime later on such simulations, since flows in bounded domains are the most natural systems. Also, you guys can try to implement the technique in the 2D code that all of us together wrote few weeks ago! You can also search independently about such methods. One helpful google-search key may be - "Volume penalization method".

One other option to simulate boundaries / bounded flows, within pseudo-spectral scheme is - instead of using the Fourier transform, one can use sine-Fourier or cosine-Fourier transforms while taking the derivatives. Such libraries are easily available within the FFTW architecture but I had a hard time in the implementation process. May be, young blood can rejuvenate the endeavor!!!

## Resources:

* [http://dx.doi.org/10.1016/j.jcp.2014.05.038](http://dx.doi.org/10.1016/j.jcp.2014.05.038)
* [https://doi.org/10.1016/j.cpc.2010.05.019](https://doi.org/10.1016/j.cpc.2010.05.019)
* [https://doi.org/10.1016/j.compfluid.2004.09.006](https://doi.org/10.1016/j.compfluid.2004.09.006)

You can get some funny videos on turbulence (some of which you can now simulate now using the 2dfft.f95 code) here.
[https://www.youtube.com/watch?v=5zI9sG3pjVU](https://www.youtube.com/watch?v=5zI9sG3pjVU)

And here is one more link [http://www.lcs-fast.com/fifth_order/](http://www.lcs-fast.com/fifth_order/) that shows some spectral simulation of air flows around a racing-car. But note that these are not pseudo-spectral. The basic idea is same though. The only difference is, instead of Fourier basis, they have used Legendre or Chebyshev basis (and sometimes Hermite also!) for better accuracy! If you are feeling confused, revise our first discussion and specifically look at the first slide!

