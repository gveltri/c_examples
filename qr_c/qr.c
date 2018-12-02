/*
  @file qr.c
  @author Gerardo Veltri
  QR Decomposition
*/
#include <stdio.h>
#include <matrix.h>
#include <stdlib.h>

int main() {
  Matrix matrix1, matrix2, matrix3;

  matrix1 = makeMatrix(9, 9, 'L', 1);
  matrix2 = makeMatrix(9, 9, 'U', 1);
  matrix3 = multiplyMatrices(matrix1, matrix2);

  draw2DMatrix(matrix1);
  draw2DMatrix(matrix2);
  draw2DMatrix(matrix3);

  freeMatrix(matrix1);
  freeMatrix(matrix2);
  freeMatrix(matrix3);

  float x = 2;

  float vec1[5] = {x,x,x,x,x};
  float vec2[5] = {x,x,x,x,x};

  x = dotProduct(5, vec1, vec2);

  printf("%.6f", x);
  return 0;
}

