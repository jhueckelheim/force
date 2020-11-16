#include <iostream>
#include <math.h>
#include "force.hpp"

template <typename T>
static T**matalloc(int n){
  auto buff = new T[n*n];
  auto a = new T*[n];
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
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      C[i][j] = 0.0;

  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < n; j++) {
	C[i][j] = C[i][j] + A[i][k] * B[k][j];
      }
    }
  }
}

// adds two nxn matrices.  C is "out" variable.
template <typename T>
static void add(int n, T**C, T**A, T**B) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = A[i][j] + B[i][j];
}

// subtracts two nxn matrices.  C is "out" variable.
template <typename T>
static void subtract(int n, T**C, T**A, T**B) {
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
static void strassenR(int n, T**C, T**A, T**B) {
  if (n <= LEAF_SIZE) {
    matrixMultiply(n, C, A, B);
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
    strassenR<T>(newSize, p1, aResult, bResult);
    // p1 = (a11+a22) * (b11+b22)

    add<T>(newSize, aResult, a21, a22); // a21 + a22
    strassenR<T>(newSize, p2, aResult, b11); // p2 = (a21+a22) * (b11)

    subtract<T>(newSize, bResult, b12, b22); // b12 - b22
    strassenR<T>(newSize, p3, a11, bResult);
    // p3 = (a11) * (b12 - b22)

    subtract<T>(newSize, bResult, b21, b11); // b21 - b11
    strassenR<T>(newSize, p4, a22, bResult);
    // p4 = (a22) * (b21 - b11)

    add<T>(newSize, aResult, a11, a12); // a11 + a12
    strassenR<T>(newSize, p5, aResult, b22);
    // p5 = (a11+a12) * (b22)

    subtract<T>(newSize, aResult, a21, a11); // a21 - a11
    add<T>(newSize, bResult, b11, b12); // b11 + b12
    strassenR<T>(newSize, p6, aResult, bResult);
    // p6 = (a21-a11) * (b11+b12)

    subtract<T>(newSize, aResult, a12, a22); // a12 - a22
    add<T>(newSize, bResult, b21, b22); // b21 + b22
    strassenR<T>(newSize, p7, aResult, bResult);
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

template <typename T>
static void fromDouble(T** matout, double** matin, int size) {
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++) {
      matout[i][j] = matin[i][j];
    }
  }
}

int main() {
  real1 **A1 = matalloc<real1>(N);
  real1 **B1 = matalloc<real1>(N);
  real1 **R1 = matalloc<real1>(N);
  real2 **A2 = matalloc<real2>(N);
  real2 **B2 = matalloc<real2>(N);
  real2 **R2 = matalloc<real2>(N);
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      A1[i][j] = (i+j);
      B1[i][j] = 1.0/(1+i+j);
      A2[i][j] = (i+j);
      B2[i][j] = 1.0/(1+i+j);
    }
  }
  matrixMultiply<real1>(N, R1, A1, B1);
  strassenR<real2>(N, R2, A2, B2);
 
  double maxrelerr = 0.0; 
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      double r1val = R1[i][j];
      double r2val = R2[i][j];
      double relerr = fabs((r1val-r2val)/r1val);
      if(relerr > maxrelerr) maxrelerr = relerr;
    }
  }
  std::cout<<"max relative err: "<<maxrelerr<<std::endl;
  matfree<real1>(N,A1);
  matfree<real1>(N,B1);
  matfree<real1>(N,R1);
  matfree<real2>(N,A2);
  matfree<real2>(N,B2);
  matfree<real2>(N,R2);
}

