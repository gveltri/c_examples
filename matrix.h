/*
  @file matrix.h
  @author Gerardo Veltri
  Base matrix functions
*/
#ifndef HEADER_FILE
#define HEADER_FILE

typedef struct _Matrix_ {
  int n, m;
  float **values;
} *Matrix;

Matrix makeMatrix(int n, int m, char type, float scalar);
void freeMatrix(Matrix matrix);
float matrixMax(Matrix matrix);
void draw2DMatrix(Matrix matrix);
Matrix multiplyMatrices(Matrix matrix1, Matrix matrix2);

#endif
