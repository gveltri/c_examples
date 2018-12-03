/*
  @file
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

Matrix copyMatrix(Matrix matrix) {
  Matrix new = malloc(sizeof(struct _Matrix_));

  new->n = matrix->n;
  new->m = matrix->m;
  new->values = malloc(new->n*sizeof(float*));
  for (int i=0;i<new->n;i++) {
    new->values[i] = malloc(new->m*sizeof(float));
    for (int j=0;j<new->m;j++) {
      new->values[i][j] = matrix->values[i][j];
    }
  }

  return new;
}

void scaleMatrix(Matrix matrix, float scalar) {
  for (int i=0;i<matrix->n;i++) {
    matrix->values[i] = malloc(matrix->m*sizeof(float));
    for (int j=0;j<matrix->m;j++) {
      matrix->values[i][j] = matrix->values[i][j] * scalar;
    }
  }
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

/* dot product */
float dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2) {
  float x = 0.0;

  switch (orient) {
  case 'R':
    for (int i=0; i<matrix1->n; i++) {
      x = x + (matrix1->values[idx1][i] * matrix2->values[idx2][i]);
    }
  case 'C':
    for (int i=0; i<matrix1->n; i++) {
      x = x + (matrix1->values[i][idx1] * matrix2->values[i][idx2]);
    }
  }

  return x;
}

float dotProductV(Matrix matrix1, Matrix matrix2) {
  float x = dotProduct('C', matrix1, 0, matrix2, 0);
  return x;
}

/* norm */
float norm(char orient, int idx, Matrix matrix) {
  float x = sqrt(dotProduct(orient, matrix, idx, matrix, idx));
  return x;
}

float normV(Matrix matrix) {
  float x = norm('C', 0, matrix);
  return x;
}

/* projection */
Matrix project(Matrix source, int idx1, Matrix target, int idx2) {
  float x;
  Matrix matrix = copyMatrix(target);
  x = (dotProduct('C', source, idx1, target, idx2) /
       dotProduct('C', target, idx1, target, idx2));
  scaleMatrix(matrix, x);
  return matrix;
}


/* matrix multiplication */
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
