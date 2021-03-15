# Code Description

## fluid_2d.f95: 

This code evaluates the two dimensional incompressible Navier-Stokes equation in vorticity-stream-function formalism using psedo-spectral method as spatial discretization and Adams-Bashforth algorithm for temporal updates. The initial condition is taken as counter-streaming flows thereby akin to Kelvin-Helmholtz type instability.

```fortran
program fft2

implicit none

integer (kind=4), parameter :: Nx = 64
integer (kind=4), parameter :: Ny = 64
integer (kind=4), parameter :: Nh = (Nx/2) + 1
real (kind=8), parameter :: pi = 3.1428

include "fftw3.f"

integer (kind=4) i,j,t
real (kind = 8) Lx, Ly, dx, dy, kx, ky, nu
real (kind = 8) time, time_min, time_max, dt
real (kind=8) x(Nx), y(Ny), omega(Nx,Ny), psi(Nx,Ny), ux(Nx,Ny), uy (Nx,Ny)
complex (kind=8) omegak(Nh,Ny), omegak_new(Nh, Ny), psik(Nh,Ny), ukx(Nh,Ny), uky(Nh,Ny)
complex (kind=8) omegak_dum(Nh,Ny),psik_dum(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny)
complex (kind=8) dt_omegak_new(Nh, Ny),  dt_omegak_old(Nh, Ny)

integer (kind=8) plan_forward, plan_backward

!====================== USER INPUTS ====================================

Lx = 2.0d0*pi
Ly = 2.0d0*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

nu = 0.10d0

time_min = 0.0d0
time_max = 10.0d0
dt = 0.0010d0

do i = 1,Nx
  x(i) = (i-1)*dx
  do j = 1,Ny
    y(j) = (j-1)*dy
    omega(i,j) = 2.0d0*dsin(x(i))*dcos(y(j))
  enddo
enddo

!=================== INITIAL TIME DATA ================================

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif  
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
    psi(i,j) = psi(i,j)/(dfloat(Nx)*dfloat(Ny))
    omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    write(10,*) x(i),y(j), omega(i,j), psi(i,j), ux(i,j), uy(i,j)
  enddo
enddo

 close(10)

!==================== MAIN PROGRAM ========================================

do time = time_min, time_max, dt

t = nint(time/dt) - int(time_min/dt)

  call derive (Nx, Ny, Nh, nu, time, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
  call ab (Nx, Ny, Nh, nu, time, dt, omegak, omegak_new, dt_omegak_old, dt_omegak_new)

do i = 1, Nh
  do j = 1, Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
    psi(i,j) = psi(i,j)/(dfloat(Nx)*dfloat(Ny))
    omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    if (mod(dfloat(t),100.0) == 0.0) then
    write(t+100,*) x(i),y(j), omega(i,j), psi(i,j), ux(i,j), uy(i,j)
    endif
  enddo
enddo
   
  if (mod(dfloat(t),100.0) == 0.0) then
  close(t+100)
  endif

enddo !time

contains

!=================== SUBROUTINE DERIVE ===================================

subroutine derive (Nx, Ny, Nh, nu, time, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
implicit none
integer (kind = 4) Nx, Ny, Nh
real (kind = 8) time, kx, ky, nu, Lx, Ly
real (kind = 8) ux(Nx,Ny), uy (Nx,Ny)
real (kind = 8) domega_dx(Nx,Ny), domega_dy(Nx,Ny) 
real (kind = 8) ux_domega_dx(Nx,Ny),uy_domega_dy(Nx,Ny)
complex (kind=8) NLkx(Nh,Ny),NLky(Nh,Ny), NLk(Nh,Ny)
complex (kind=8) omegak(Nh,Ny), omegak_new(Nh,Ny), psik(Nh,Ny), ukx(Nh,Ny), uky(Nh,Ny)
complex (kind=8) i_kx_omegak(Nh,Ny),i_ky_omegak(Nh,Ny)
complex (kind=8) omegak_dum(Nh,Ny),psik_dum(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny)
complex (kind=8) dt_omegak_new(Nh, Ny),  dt_omegak_old(Nh, Ny)

Lx = 2.0d0*pi
Ly = 2.0d0*pi

do i = 1, Nh
  do j = 1, Ny
  omegak_dum(i,j) = omegak(i,j)
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0) 
      else
      psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0) * kx * omegak(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0) * ky * omegak(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0) * kx * omegak(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0) * ky * omegak(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak, domega_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1,Nx
  do j = 1,Ny
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dx(i,j) = domega_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dy(i,j) = domega_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_domega_dx(i,j) = ux(i,j) * domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j) * domega_dy(i,j)
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! De-Aliazing Technique - 2/3 Truncation...
do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    if ( dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 ) then
    NLkx(i,j) = 0.0d0
    NLky(i,j) = 0.0d0
    endif
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    if ( dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 ) then
    NLkx(i,j) = 0.0d0
    NLky(i,j) = 0.0d0
    endif
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk(i,j) = NLkx(i,j) + NLky(i,j)
    NLk(i,j) = NLk(i,j) + nu * (kx*kx + ky*ky) * omegak(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    NLk(i,j) = NLkx(i,j) + NLky(i,j)
    NLk(i,j) = NLk(i,j) + nu * (kx*kx + ky*ky) * omegak(i,j)
  enddo
enddo

do i = 1,Nh
  do j = 1,Ny
    dt_omegak_new(i,j) = -NLk(i,j)
  enddo
enddo

return

end subroutine derive

!=================== SUBROUTINE ADAMS-BASHFORTH ==========================

subroutine ab (Nx, Ny, Nh, nu, time, dt, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
implicit none
integer (kind = 4) Nx, Ny, Nh
real (kind = 8) time, dt, nu
complex (kind=8) omegak(Nh,Ny), omegak_new(Nh, Ny)
complex (kind=8) dt_omegak_new(Nh, Ny),  dt_omegak_old(Nh, Ny)

do i = 1,Nh
  do j = 1, Ny
    omegak_new(i,j) = omegak(Nh,Ny) + ( (3.0d0/2.0d0) * dt_omegak_new(i,j) - (1.0d0/2.0d0) * dt_omegak_old(i,j)) * dt
  enddo
enddo

end subroutine ab

!========================================================================
end program fft2
```
