/*
  @file matrix.c
  @author Gerardo Veltri
  Base matrix functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matrix.h>

/* constants for rendering tables */
const int PADDING = 1;
const int PRECISION = 2;
const char *FORMATTING  = "%.2f";


void allocValues(float value, char type, Matrix matrix) {
  for (int i=0;i<matrix->n;i++) {
    matrix->values[i] = malloc(matrix->m*sizeof(float));
    for (int j=0;j<matrix->m;j++) {

      /* base matrix types */
      switch (type) {
      case 'V':
        matrix->values[i][j] = value;
        break;
      case 'I':
        if (i == j)
          matrix->values[i][j] = value;
        else
          matrix->values[i][j] = 0;
        break;
      case 'U':
        if (i <= j)
          matrix->values[i][j] = value;
        else
          matrix->values[i][j] = 0;
        break;
      case 'L':
        if (i >= j)
          matrix->values[i][j] = value;
        else
          matrix->values[i][j] = 0;
        break;
      }

    }
  }
}

Matrix makeMatrix(int n, int m, char type, float scalar) {
  Matrix matrix = malloc(sizeof(struct _Matrix_));

  matrix->n = n;
  matrix->m = m;
  matrix->values = malloc(n*sizeof(float*));
  allocValues(scalar, type, matrix);

  return matrix;
}

void freeMatrix(Matrix matrix) {
  for (int i; i<matrix->n; i++) {
    free(matrix->values[i]);
  }

  free(matrix->values);
  free(matrix);
}

float maxValue(int m, float *array) {
  int max = array[0];

  for (int i=1; i < m; i++) {
    if (max < array[i])
      max = array[i];
  }

  return max;
}

float dotProduct(int n, float *vector1, float *vector2) {
  float x = 0;

  for (int i=0; i<n; i++) {
    x = x + (vector1[i] * vector2[i]);
  }


  return x;
}

float matrixMax(Matrix matrix) {
  float max = maxValue(matrix->m, matrix->values[0]);

  for (int i=1; i < matrix->n; i++) {
    float curr = maxValue(matrix->m, matrix->values[i]);
    if (max < curr)
      max = curr;
  }

  return max;
}

int numDigits(int n) {
  int add = 0;
  if (n < 0)
    n = n * -1;
    add = 1;

  if (n < 1)
    n = 1;
  else
    n = floor(log10(n)) + 1;

  return n;
}

int numDigitsFloat(float x) {
  int n = floor(x);

  n = numDigits(n);

  return n;
}

void draw2DMatrix(Matrix matrix) {
  float matrix_max = matrixMax(matrix);
  int max_width = numDigitsFloat(matrix_max);
  max_width = max_width + PADDING + PRECISION - 1;

   /* header */
  printf("\n");
  printf("(%d, %d) ", matrix->n, matrix->m);
  int offset = numDigits(matrix->n) + numDigits(matrix->m) + 5;

  for (int i=offset; i<(max_width*matrix->m); i++) {
    printf("-");
  }

  for (int i=0; i<matrix->n; i++) {
    printf("\n");
    for (int j=0; j<PADDING; j++) {
      printf("\n");
    }
    for (int j=0; j<matrix->m; j++) {
      float value = matrix->values[i][j];
      printf(FORMATTING, value);
      int blank_fill = max_width - numDigitsFloat(value);
      for (int i=0; i<blank_fill; i++) {
        printf(" ");
      }
    }
  }

  printf("\n");
  for (int i=0; i<(max_width*matrix->m); i++) {
    printf("-");
  }
  for (int i=0; i<3; i++) printf("\n");
}

Matrix multiplyMatrices(Matrix matrix1, Matrix matrix2) {
  Matrix matrix3;

  matrix3 = makeMatrix(matrix1->n, matrix2->m, 'V', 0);

  for (int i=0; i<matrix1->n; i++) {
    for (int j=0; j<matrix2->m; j++) {
      float value = 0;
      for (int k=0; k<matrix1->m; k++) {
        value = value + (matrix1->values[i][k] * matrix2->values[k][j]);
      }
      matrix3->values[i][j] = value;
    }
  }
  return matrix3;
}
