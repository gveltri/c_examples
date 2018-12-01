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
  matrix1 = makeMatrix(5,5,'V', 1);
  matrix2 = makeMatrix(5,1,'I', 1);
  matrix3 = multiplyMatrices(matrix1, matrix2);
  draw2DMatrix(matrix1);
  draw2DMatrix(matrix2);
  draw2DMatrix(matrix3);
  freeMatrix(matrix1);
  freeMatrix(matrix2);
  freeMatrix(matrix3);
  return 0;
}

