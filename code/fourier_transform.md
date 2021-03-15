```fortran
PROGRAM Fourier_Transform

implicit none

real*8, parameter :: pi = 3.14159265358979323846d0

integer, parameter :: N = 1000 ! Number of data points.

INTEGER, PARAMETER :: UN = N/2+1 ! You may use lower mode-numbers to check convergence.

real*8, parameter :: Length = 6.0d0*pi ! Maximum value of Independent variable.

real*8 dx,f,x0,a0,a,b,k,ft,Re,Im

real*8 x(N),y(0:N),y_rec(0:N),rfft(N),ifft(N)

integer i,u

external f

!open(unit=1,file='Data_File.dat',status='old')

open(unit=10,file='Input_Function.dat',status='unknown')

open(unit=20,file='Fourier_Coeff.dat',status='unknown')

open(unit=30,file='Input_Function_Recovered.dat',status='unknown')

x0 = 0.0d0

dx = Length/dfloat(N)

y(0) = f(x0)

do i = 1,N-1

  x(i) = x0 + dx*dfloat(i)

  y(i) = f(x(i))

!  read(1,*) x(i), y(i)

  write(10,*) x(i), y(i)

enddo ! i

a0 = y(0)

do i = 1,N-1

  a0 = a0 + y(i)

enddo ! i

write(20,*) 0, a0/dfloat(N)

do u = 1,UN-1

  k = 2.0d0*pi*dfloat(u)/Length

  rfft(u) = f(x0)

  ifft(u) = 0.0d0

    do i = 1,N-1

      a = 2.0d0*pi*x(i)/Length

      rfft(u) = rfft(u) + y(i)*dcos(dfloat(u)*a)

      ifft(u) = ifft(u) + y(i)*dsin(dfloat(u)*a)

    enddo ! i

  b = dsqrt(rfft(u)**2 + ifft(u)**2)/dfloat(N)

  ft = 2.0d0*b

  write(20,*) k, ft!, rfft(u)/dfloat(N), ifft(u)/dfloat(N)

enddo ! u

do u = N-1,(N-UN)+1,-1

!  k = 2.0d0*pi*dfloat(u)/Length

  rfft(u) = f(x0)

  ifft(u) = 0.0d0

    do i = 1,N-1

      a = 2.0d0*pi*x(i)/Length

      rfft(u) = rfft(u) + y(i)*dcos(dfloat(u)*a)

      ifft(u) = ifft(u) + y(i)*dsin(dfloat(u)*a)

    enddo ! i

!  b = dsqrt(rfft(u)**2 + ifft(u)**2)/dfloat(N)

!  ft = 2.0d0*b

!  write(20,*) k, ft!, rfft(u)/dfloat(N), ifft(u)/dfloat(N)

enddo ! u

do i = 1,N-1

  Re = 0.0d0

  Im = 0.0d0

    do u = 1,UN-1

      Re = Re + rfft(u)*dcos(2.0d0*pi*dfloat(u)*x(i)/Length)

      Im = Im + ifft(u)*dsin(2.0d0*pi*dfloat(u)*x(i)/Length)

    enddo ! u

    do u = N-1,(N-UN)+1,-1

      Re = Re + rfft(u)*dcos(2.0d0*pi*dfloat(u)*x(i)/Length)

      Im = Im + ifft(u)*dsin(2.0d0*pi*dfloat(u)*x(i)/Length)

    enddo

  y_rec(i) = (a0 + Re + Im)/dfloat(N)

  write(30,*) x(i), y_rec(i)

enddo ! i

END PROGRAM Fourier_Transform

!===========================================================

FUNCTION f(x)

IMPLICIT NONE

real*8, parameter :: pi = 3.14159265358979323846d0

REAL*8 f,x,lamda

!

! YOUR INPUT FUNCTION HERE

!

lamda = 2.0d0

!f= dsin(5.0d0*x)

!f = 7.0d0 + 5.0d0*dsin(2.0d0*x)+ 3.0d0*dcos(8.0d0*x/lamda)

if ((x .ge. 1.0d0) .and. (x .le. 5.0d0)) then

  f = 4.0d0

  else

  f = -1.0d0

endif

END FUNCTION
```
