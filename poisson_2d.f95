program fft2

implicit none

integer (kind=4), parameter :: Nx = 64
integer (kind=4), parameter :: Ny = 64
integer (kind=4), parameter :: Nh = (Nx/2) + 1
real (kind=8), parameter :: pi = 3.1428

include "fftw3.f"

integer (kind=4) i,j
real (kind = 8) Lx, Ly, dx, dy, kx, ky
real (kind=8) x(Nx), y(Ny), rho(Nx,Ny), phi(Nx,Ny)
complex (kind=8) rhok(Nh,Ny), phik(Nh,Ny)

integer (kind=8) plan_forward, plan_backward

Lx = 2.0d0*pi
Ly = 2.0d0*pi
dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

do i = 1,Nx
x(i) = (i-1)*dx
  do j = 1,Ny
  y(j) = (j-1)*dy
  rho(i,j) = 2.0d0*dsin(x(i))*dcos(y(j))
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho, rhok, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      phik(i,j) = (0.0d0,0.0d0)
      else
      phik(i,j) = rhok(i,j)/(kx*kx+ky*ky)
      endif  
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    phik(i,j) = rhok(i,j)/(kx*kx+ky*ky)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, phik, phi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
  write(10,*) x(i),y(j), rho(i,j), phi(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo

end program fft2
