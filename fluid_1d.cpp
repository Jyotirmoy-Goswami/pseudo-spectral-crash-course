/* Modified: burgers2.cpp */
/* Solving 1D Burger equation using using fftw library */
/*This program is written by Gunjan Sharma, CPP-IPR, India on 16 July, 2020*/

#include<iostream>
#include<fftw3.h>
#include<cmath>
#include<complex>
#define pi 3.1415

int main()
{
   int i,j;
   int N = 1024;
   int Nh = N/2 + 1;
   double x[N], u[N], u_dum[N], u_rec[N], uk_abs[Nh], du[N], force_o[N], force_n[N], ddu[N];
   fftw_complex uk[Nh], uk_dum[Nh], ik_uk[Nh], k2_uk[Nh];
   fftw_plan f;
   fftw_plan b1;
   fftw_plan b2;
   double L = 2*pi, k, nu = 0.0001;
   double dx = L/double(N);
   double dt = 0.00010;
   double v = dx/dt;
   double time=0.0, tmax=1.20;


   FILE *file1;
   file1 = fopen("bur_gauss1.txt","w");


   for(i=0;i<N;i++)
   {
      x[i] = double(i)*dx;
      u[i] = sin(x[i]); //exp(-(x[i]-pi)*(x[i]-pi));
      u_dum[i] = u[i];
   }

   /* time loop */
   while(time<tmax)
   {

   j=int(time/dt);

   f = fftw_plan_dft_r2c_1d(N, u_dum, uk, FFTW_ESTIMATE);
   fftw_execute(f);
   fftw_destroy_plan(f);
   fftw_cleanup();


   for(i=0;i<Nh;i++)
   {
      k = 2*pi*i/L;

      // find du/dx
      ik_uk[i][0] = -k*uk[i][1];
      ik_uk[i][1] = k*uk[i][0];

      // find d2u/dx2
      k2_uk[i][0] = -k*k*uk[i][0];
      k2_uk[i][1] = -k*k*uk[i][1];
   }

   b1 = fftw_plan_dft_c2r_1d(N, ik_uk, du, FFTW_ESTIMATE);
   fftw_execute(b1);
   fftw_destroy_plan(b1);
   fftw_cleanup();

   b2 = fftw_plan_dft_c2r_1d(N, k2_uk, ddu, FFTW_ESTIMATE);
   fftw_execute(b2);
   fftw_destroy_plan(b2);
   fftw_cleanup();

   for(i=0;i<N;i++)
   {
      du[i] = du[i]/double(N);
      ddu[i] = ddu[i]/double(N);

      if(j%3000==0){
         fprintf(file1, "%d %g %g\n",i,x[i],u[i]);
         std::cout << u[i] <<  std::endl;
         }
   }

   for(i=0;i<N;i++)
   {
     force_n[i] =  -u[i]*du[i] + nu*du[i] ;

      if(j==0) u[i] = u[i] + force_n[i]*dt;
      else u[i] = u[i] + dt*( (3.0/2.0)*force_n[i] - (1.0/2.0)*force_o[i] );

     force_o[i] = force_n[i];
     u_dum[i] = u[i];
   }
   time += dt;
   }


   fclose(file1);

   return 0;
}
