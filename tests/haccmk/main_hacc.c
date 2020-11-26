#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "mpfrcpp_tpl.h"
#include "force.hpp"
#ifdef CENA
typedef freal<float> usereal;
#elif CDBL
typedef freal<double> usereal;
#elif CLDB
typedef freal<long double> usereal;
#elif MPFR
typedef mpfrcpp<MPFRPR> usereal;
#elif QUAD
#include <quadmath.h>
typedef __float128 usereal;
usereal sqrt(usereal val) {
  return sqrtq(val);
}
std::ostream& operator<<(std::ostream &ost, const usereal &ad){
   ost << (double)ad;
   return ost;
}
#elif DBLE
typedef double usereal;
#elif LDBL
typedef long double usereal;
#else
typedef float usereal;
#endif

void Step10_orig( int count1, usereal xxi, usereal yyi, usereal zzi, usereal fsrrmax2, usereal mp_rsm2, usereal *xx1, usereal *yy1, usereal *zz1, usereal *mass1, usereal *dxi, usereal *dyi, usereal *dzi )
{

    const usereal ma0 = 0.269327, ma1 = -0.0750978, ma2 = 0.0114808, ma3 = -0.00109313, ma4 = 0.0000605491, ma5 = -0.00000147177;
    
    usereal dxc, dyc, dzc, m, r2, f, xi, yi, zi;
    int j;

    xi = 0.; yi = 0.; zi = 0.;
    usereal zero = 0.0, one = 1.0;

    for ( j = 0; j < count1; j++ ) 
    {
        dxc = xx1[j] - xxi;
        dyc = yy1[j] - yyi;
        dzc = zz1[j] - zzi;
  
        r2 = dxc * dxc + dyc * dyc + dzc * dzc;
       
        m = ( r2 < fsrrmax2 ) ? mass1[j] : zero;

        f = (r2 + mp_rsm2);
        f = one/(f*sqrt(f));

        f =  f - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        
        f = ( r2 > zero ) ? m * f : zero;

        xi = xi + f * dxc;
        yi = yi + f * dyc;
        zi = zi + f * dzc;
    }

    *dxi = xi;
    *dyi = yi;
    *dzi = zi;
}

#ifdef TIMEBASE
extern unsigned long long timebase();
#else
extern double mysecond();
#endif

// Frequency in MHz, L1 cache size in bytes, Peak MFlops per node */
// BG/Q 
#define MHz 1600e6
#define NC 16777216
#define PEAK (16.*12800.)

// 4 core Intel(R) Xeon(R) CPU E5-2643 0 @ 3.30GHz - Rick's machine celero 
// Peak flop rate: AVX: 8 Flops/cycle * 4 cores * 3291.838 Mcycles/second 
/*
#define MHz 3291.838e6     
#define NC (64*1024*1024)
#define PEAK 105338.816
*/

// 4 core Intel(R) Xeon(R) CPU           E5430  @ 2.66GHz - crush
// Peak flop rate: SSE: 4 Flops/cycle * 4 cores * 2666.666 Mcycles/second
/*
#define MHz 2666.7e6
#define NC ( 256*1024*1024 )
#define PEAK 42667.2
*/

// 4 core Intel(R) Core(TM) i7-3820 CPU @ 3.60GHz - strength
// Peak flop rate: AVX: 8 Flops/cycle * 4 cores * 3600 MCycles/second
/*
#define MHz 3600.e6
#define NC ( 32*1024*1024 )
#define PEAK 115200.
*/

#define N 15000      /* Vector length, must be divisible by 4  15000 */

#define ETOL 1.e-4  /* Tolerance for correctness */

int main( int argc, char **argv )
{
  static usereal xx[N], yy[N], zz[N], mass[N], vx1[N], vy1[N], vz1[N];
  usereal fsrrmax2, mp_rsm2, fcoeff, dx1, dy1, dz1;

  char  M1[NC], M2[NC];
  int n, count, i, rank, nprocs;
  double tm1, tm2, tm3, tm4, total = 0;
  double t3, elapsed = 0.0;
  usereal final, validation;

  //MPI_Init( &argc, &argv );
  //MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  //MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  
  rank = 0;
  nprocs = 1;

  count = 327;

  if ( rank == 0 ) 
  {  
     printf( "count is set %d\n", count );
     printf( "Total MPI ranks %d\n", nprocs );
  } 

#pragma omp parallel
{
  if ( (rank == 0) && (omp_get_thread_num() == 0) )
  {
     printf( "Number of OMP threads %d\n\n", omp_get_num_threads() );
     //printf( "      N         Time,us        Validation result\n" );
  }   
}

#ifdef TIMEBASE
  tm3 = timebase();
#endif

  final = 0.;

  for ( n = 400; n < N; n = n + 20 ) 
  {
      /* Initial data preparation */
      fcoeff = 0.23f;  
      fsrrmax2 = 0.5f; 
      mp_rsm2 = 0.03f;
      dx1 = 1.0f/(float)n;
      dy1 = 2.0f/(float)n;
      dz1 = 3.0f/(float)n;
      xx[0] = 0.f;
      yy[0] = 0.f;
      zz[0] = 0.f;
      mass[0] = 2.f;
      
      for ( i = 1; i < n; i++ )
      {
          usereal ir = 0.01*i;
          xx[i] = xx[i-1] + dx1;
          yy[i] = yy[i-1] + dy1;
          zz[i] = zz[i-1] + dz1;
          mass[i] = ir + xx[i];
      }
    
      for ( i = 0; i < n; i++ )
      {
          vx1[i] = 0.f;
          vy1[i] = 0.f;
          vz1[i] = 0.f;
      }
    
      /* Data preparation done */
    
    
      /* Clean L1 cache */
      for ( i = 0; i < NC; i++ ) M1[i] = 4;
      for ( i = 0; i < NC; i++ ) M2[i] = M1[i];

#ifdef TIMEBASE
      tm1 = timebase();
#else
      tm1 = mysecond();
#endif

      #pragma omp parallel for private( dx1, dy1, dz1 )
      for ( i = 0; i < count; ++i)
      {
        Step10_orig( n, xx[i], yy[i], zz[i], fsrrmax2, mp_rsm2, xx, yy, zz, mass, &dx1, &dy1, &dz1 );
    
        vx1[i] = vx1[i] + dx1 * fcoeff;
        vy1[i] = vy1[i] + dy1 * fcoeff;
        vz1[i] = vz1[i] + dz1 * fcoeff;
      }

#ifdef TIMEBASE
      tm2 = timebase();
#else
      tm2 = mysecond();
#endif

      validation = 0.;
      for ( i = 0; i < n; i++ )
      {
         validation = validation + ( vx1[i] + vy1[i] + vz1[i] );
      }

      final = final + validation;
    
#ifdef TIMEBASE      
      t3 = 1e6 * (double)(tm2 - tm1) / MHz; // time in us
#else
      t3 = (tm2 - tm1) * 1e6;
#endif 

      elapsed = elapsed + t3;
    
  }

#ifdef TIMEBASE
  tm4 = timebase();
  if ( rank == 0 )
  {
      printf( "\nKernel elapsed time, s: %18.8lf\n", elapsed*1e-6 );
      printf(   "Total  elapsed time, s: %18.8lf\n", (double)(tm4 - tm3) / MHz ); 
      printf(   "Result validation: %18.8lf\n", final );
      printf(   "Result expected  : 6636045675.12190628\n" );
  }
#else
  if ( rank == 0 )
  {
      printf( "\nKernel elapsed time, s: %18.8lf\n", elapsed*1e-6 );
      printf(   "Result validation: " );
      std::cout<<std::setprecision(36)<<final<<std::endl;
      std::cout<<"Result expected:   1.03733e+07"<<std::endl;
      //printf(   "Result expected  : 6636045675.12190628\n" );
  }    
#endif
  
  //MPI_Finalize();

  return 0;
}

