# Quizzes for Lecture - 6

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

[Back to Lecture - 6](https://github.com/RupakMukherjee/fluid_teaching#lecture-6)
