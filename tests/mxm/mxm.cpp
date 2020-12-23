#include <iostream>
#include <math.h>
#include <chrono>
#include "force.hpp"
#include "randomhelper.hpp"
#include "mpfrcpp_tpl.h"

typedef mpfrcpp<200> realref;

template <typename T>
static T**matalloc(int n){
  auto buff = new T[n*n];
  auto a = new T*[n];
  #pragma omp parallel for
  for(int i = 0; i < n; ++i)
    a[i] = &buff[n*i];
  return a;
}

template <typename T>
static void matfree(int n, T** mat) {
  delete[] mat[0];
  delete [] mat;
}

// naive matrix multiplication, no Strassen here.
template <typename T>
static void matrixMultiply(int n, T**C, T**A, T**B) {
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
	C[i][j] = C[i][j] + A[i][k] * B[k][j];
      }
    }
  }
}

// adds two nxn matrices.  C is "out" variable.
template <typename T>
static void add(int n, T**C, T**A, T**B) {
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = A[i][j] + B[i][j];
}

// subtracts two nxn matrices.  C is "out" variable.
template <typename T>
static void subtract(int n, T**C, T**A, T**B) {
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = A[i][j] - B[i][j];
}

// Strassen algorithm from
// https://martin-thoma.com/strassen-algorithm-in-python-java-cpp/
// I'm just going to assume n is a power of 2.
// There is no problem dealing with the general case but need more
// time!

// multiplies two nxn matrices, storing result in C
template <typename T>
static void strassenR(int n, T**C, T**A, T**B, int leafsize) {
  if (n <= leafsize) {
    matrixMultiply<T>(n, C, A, B);
  } else {
    // initializing the new sub-matrices
    int newSize = n / 2;
    T**a11 = matalloc<T>(newSize);
    T**a12 = matalloc<T>(newSize);
    T**a21 = matalloc<T>(newSize);
    T**a22 = matalloc<T>(newSize);

    T**b11 = matalloc<T>(newSize);
    T**b12 = matalloc<T>(newSize);
    T**b21 = matalloc<T>(newSize);
    T**b22 = matalloc<T>(newSize);

    T**aResult = matalloc<T>(newSize);
    T**bResult = matalloc<T>(newSize);

    T**p1 = matalloc<T>(newSize);
    T**p2 = matalloc<T>(newSize);
    T**p3 = matalloc<T>(newSize);
    T**p4 = matalloc<T>(newSize);
    T**p5 = matalloc<T>(newSize);
    T**p6 = matalloc<T>(newSize);
    T**p7 = matalloc<T>(newSize);

    T**c11 = matalloc<T>(newSize);
    T**c12 = matalloc<T>(newSize);
    T**c21 = matalloc<T>(newSize);
    T**c22 = matalloc<T>(newSize);

    // dividing the matrices in 4 sub-matrices:
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < newSize; i++) {
      for (int j = 0; j < newSize; j++) {
	a11[i][j] = A[i][j]; // top left
	a12[i][j] = A[i][j + newSize]; // top right
	a21[i][j] = A[i + newSize][j]; // bottom left
	a22[i][j] = A[i + newSize][j + newSize]; // bottom right
	
	b11[i][j] = B[i][j]; // top left
	b12[i][j] = B[i][j + newSize]; // top right
	b21[i][j] = B[i + newSize][j]; // bottom left
	b22[i][j] = B[i + newSize][j + newSize]; // bottom right
      }
    }
    // Calculating p1 to p7:
    add<T>(newSize, aResult, a11, a22);
    add<T>(newSize, bResult, b11, b22);
    strassenR<T>(newSize, p1, aResult, bResult, leafsize);
    // p1 = (a11+a22) * (b11+b22)

    add<T>(newSize, aResult, a21, a22); // a21 + a22
    strassenR<T>(newSize, p2, aResult, b11, leafsize); // p2 = (a21+a22) * (b11)

    subtract<T>(newSize, bResult, b12, b22); // b12 - b22
    strassenR<T>(newSize, p3, a11, bResult, leafsize);
    // p3 = (a11) * (b12 - b22)

    subtract<T>(newSize, bResult, b21, b11); // b21 - b11
    strassenR<T>(newSize, p4, a22, bResult, leafsize);
    // p4 = (a22) * (b21 - b11)

    add<T>(newSize, aResult, a11, a12); // a11 + a12
    strassenR<T>(newSize, p5, aResult, b22, leafsize);
    // p5 = (a11+a12) * (b22)

    subtract<T>(newSize, aResult, a21, a11); // a21 - a11
    add<T>(newSize, bResult, b11, b12); // b11 + b12
    strassenR<T>(newSize, p6, aResult, bResult, leafsize);
    // p6 = (a21-a11) * (b11+b12)

    subtract<T>(newSize, aResult, a12, a22); // a12 - a22
    add<T>(newSize, bResult, b21, b22); // b21 + b22
    strassenR<T>(newSize, p7, aResult, bResult, leafsize);
    // p7 = (a12-a22) * (b21+b22)

    // calculating c21, c21, c11 e c22:
    add<T>(newSize, c12, p3, p5); // c12 = p3 + p5
    add<T>(newSize, c21, p2, p4); // c21 = p2 + p4

    add<T>(newSize, aResult, p1, p4); // p1 + p4
    add<T>(newSize, bResult, aResult, p7); // p1 + p4 + p7
    subtract<T>(newSize, c11, bResult, p5);
    // c11 = p1 + p4 - p5 + p7

    add<T>(newSize, aResult, p1, p3); // p1 + p3
    add<T>(newSize, bResult, aResult, p6); // p1 + p3 + p6
    subtract<T>(newSize, c22, bResult, p2);
    // c22 = p1 + p3 - p2 + p6

    // Grouping the results obtained in a single matrix:
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < newSize; i++) {
      for (int j = 0; j < newSize; j++) {
	C[i][j] = c11[i][j];
	C[i][j + newSize] = c12[i][j];
	C[i + newSize][j] = c21[i][j];
	C[i + newSize][j + newSize] = c22[i][j];
      }
    }

    matfree<T>(newSize,a11);
    matfree<T>(newSize,a12);
    matfree<T>(newSize,a21);
    matfree<T>(newSize,a22);

    matfree<T>(newSize,b11);
    matfree<T>(newSize,b12);
    matfree<T>(newSize,b21);
    matfree<T>(newSize,b22);

    matfree<T>(newSize,aResult);
    matfree<T>(newSize,bResult);

    matfree<T>(newSize,p1);
    matfree<T>(newSize,p2);
    matfree<T>(newSize,p3);
    matfree<T>(newSize,p4);
    matfree<T>(newSize,p5);
    matfree<T>(newSize,p6);
    matfree<T>(newSize,p7);

    matfree<T>(newSize,c11);
    matfree<T>(newSize,c12);
    matfree<T>(newSize,c21);
    matfree<T>(newSize,c22);
  }
}

template <typename T, typename Tref, typename Tin>
void benchmark(int N, int leafsize, Tin** Ain, Tin** Bin, Tref** Rref) {
  T **A = matalloc<T>(N);
  T **B = matalloc<T>(N);
  T **R = matalloc<T>(N);
  {
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        A[i][j] = Ain[i][j];
        B[i][j] = Bin[i][j];
      }
    }
    auto t_start = std::chrono::high_resolution_clock::now();
    matrixMultiply<T>(N, R, A, B);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    Tref relerrsum = 0.0; 
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        Tref val = (Tref)R[i][j];
        Tref refval = Rref[i][j];
        Tref relerr = (val-refval)/refval;
        relerrsum = relerrsum + fabs(relerr);
      }
    }
    std::cout<<"t_mxm "<<time<<" err "<<relerrsum<<std::endl;
  }
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    strassenR<T>(N, R, A, B, leafsize);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    Tref relerrsum = 0.0; 
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        Tref val = (Tref)R[i][j];
        Tref refval = Rref[i][j];
        Tref relerr = (val-refval)/refval;
        relerrsum = relerrsum + fabs(relerr);
      }
    }
    std::cout<<"         t_str "<<time<<" err "<<relerrsum<<std::endl;
  }
  matfree<T>(N,A);
  matfree<T>(N,B);
  matfree<T>(N,R);
}

int main(int argc, char** argv) {
  if(argc < 2) {
    std::cout<<"Usage: ./exec <N> [LEAFSIZE]"<<std::endl;
    exit(-1);
  }
  int sleafsize, dleafsize, ldleafsize, qleafsize, csleafsize, cdleafsize, cldleafsize;
  if(argc == 3) {
    sleafsize = atoi(argv[2]);
    dleafsize = sleafsize;
    ldleafsize = sleafsize;
    qleafsize = sleafsize;
    csleafsize = sleafsize;
    cdleafsize = sleafsize;
    cldleafsize = sleafsize;
  }
  else {
    sleafsize = 256;
    dleafsize = 256;
    ldleafsize = 512;
    qleafsize = 128;
    csleafsize = 256;
    cdleafsize = 128;
    cldleafsize = 512;
  }
  int N = atoi(argv[1]);
  std::cout<<argv[0]<<" N "<<N<<std::endl;

  realref **Aref = matalloc<realref>(N);
  realref **Bref = matalloc<realref>(N);
  realref **Rref = matalloc<realref>(N);
  __float128 **Aq = matalloc<__float128>(N);
  __float128 **Bq = matalloc<__float128>(N);
  __float128 **Rq = matalloc<__float128>(N);
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      Aq[i][j] = rand_quad();
      Bq[i][j] = rand_quad();
      Aref[i][j] = (realref)Aq[i][j];
      Bref[i][j] = (realref)Bq[i][j];
    }
  }
  strassenR<realref>(N, Rref, Aref, Bref, 128);

  std::cout<<"single   "; benchmark<float, realref, __float128>(N, sleafsize, Aq, Bq, Rref);
  std::cout<<"double   "; benchmark<double, realref, __float128>(N, dleafsize, Aq, Bq, Rref);
  std::cout<<"ldouble  "; benchmark<long double, realref, __float128>(N, ldleafsize, Aq, Bq, Rref);
  std::cout<<"quad     "; benchmark<__float128, realref, __float128>(N, qleafsize, Aq, Bq, Rref);
  std::cout<<"csingle  "; benchmark<freal<float>, realref, __float128>(N, csleafsize, Aq, Bq, Rref);
  std::cout<<"cdouble  "; benchmark<freal<double>, realref, __float128>(N, cdleafsize, Aq, Bq, Rref);
  std::cout<<"cldouble "; benchmark<freal<long double>, realref, __float128>(N, cldleafsize, Aq, Bq, Rref);

  matfree<realref>(N,Aref);
  matfree<realref>(N,Bref);
  matfree<realref>(N,Rref);
  matfree<__float128>(N,Aq);
  matfree<__float128>(N,Bq);
  matfree<__float128>(N,Rq);
}

