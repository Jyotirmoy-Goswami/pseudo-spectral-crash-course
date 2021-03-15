```fortran
program fft2

implicit none

integer (kind=4), parameter :: Nx = 64
integer (kind=4), parameter :: Ny = 64
integer (kind=4), parameter :: Nh = (Nx/2) + 1
real (kind=8), parameter :: pi = 3.1428

include "fftw3.f"

integer (kind=4) i,j
real (kind=8) ux(Nx,Ny),ux_dum(Nx,Ny)
complex (kind=8) ukx(Nh,Ny), ukx_dum(Nh,Ny)

integer (kind=8) plan_forward, plan_backward

Lx = 2.0d0*pi
Ly = 2.0d0*pi
dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

do i = i,Nx
x(i) = (i-1)*dx
  do j = 1,Ny
  y(j) = (j-1)*dy
  rho(i,j) = dsin(x(i))*dcos(y(j))
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux, ukx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

do i = 1,Nh
  do j = 1,Ny
  ukx_dum(i,j) = ukx(i,j)
  write(*,*) i,j,ukx_dum(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux_dum, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
  write(*,*) i,j, ux(i,j), ux_dum(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo

end

```
