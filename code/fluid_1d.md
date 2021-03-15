```fortran
program burgers
implicit none

include "fftw3.f"

integer, parameter :: N = 1024
integer, parameter :: Nh = N/2+1
double precision, parameter :: pi = 3.1428

double precision :: x, ux, dux, ddux
dimension x(N), ux(N), dux(N), ddux(N)
double complex :: uk,duk,dduk
dimension uk(Nh),duk(Nh),dduk(Nh)
integer*8 :: plan
real*8 :: L

integer :: i,t
real*8 dx,k,time,time_min,time_max,dt,nu
real*8, dimension (N) :: force_u_n, force_u_o

L = 2.0d0*pi
dx = L/dfloat(N)

time_min = 0.0d0
time_max = 2.0d0
dt = 0.00010d0

nu = 0.01d0

do i = 1,N
  x(i) = dfloat(i-1)*dx
  ux(i) = dsin(x(i))
enddo
 
do time = time_min, time_max, dt

t = nint(time/dt) - int(time_min/dt)

call dfftw_plan_dft_r2c_1d(plan, N, ux, uk, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, ux, uk)
call dfftw_destroy_plan(plan)

do i = 1,Nh
  k = 2.0d0*pi*dfloat(i-1)/L
  duk(i) = (0.0d0,1.0d0) * k * uk(i)
  dduk(i) = - k * k * uk(i)
enddo

call dfftw_plan_dft_c2r_1d(plan, N, duk, dux, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, duk, dux)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan, N, dduk, ddux, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, dduk, ddux)
call dfftw_destroy_plan(plan)

do i = 1,N
dux(i) = dux(i)/dfloat(N)
ddux(i) = ddux(i)/dfloat(N)
  if (t .ge. 1000 .and. mod(t,500) == 0) then
  write(t,*) t, x(i), ux(i)
  endif
enddo

do i = 1,N
  force_u_n(i) = -ux(i)*dux(i) + nu*ddux(i)
    if (t==0) then
    ux(i) = ux(i) + dt * force_u_n(i)
    else
    ux(i) = ux(i) + dt* ( (3.0/2.0) * force_u_n(i) - (1.0/2.0) * force_u_o(i))
    endif
 force_u_o(i) = force_u_n(i)
enddo

enddo

end program burgers
```
