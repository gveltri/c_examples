/*
  @file
  @author Gerardo Veltri
  QR Decomposition
*/
#include <stdio.h>
#include <matrix.h>
#include <stdlib.h>

int main() {
  Matrix matrix1, matrix2, matrix3;

  matrix1 = makeMatrix(4, 1, 'I', 1);
  matrix2 = makeMatrix(4, 1, 'V', 1);

  draw2DMatrix(matrix1);
  draw2DMatrix(matrix2);

  float x = dotProductV(matrix1, matrix2);

  printf("%.6f", x);

  printf("\n");

  float y = normV(matrix1);

  printf("%.6f", y);

  matrix3 = project(matrix1, 0, matrix2, 0);

  draw2DMatrix(matrix3);

  freeMatrix(matrix1);
  freeMatrix(matrix2);

  return 0;
}

