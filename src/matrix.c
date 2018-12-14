/*
  @file matrix.c
  @author Gerardo Veltri
  Base matrix functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mem.h>
#include <matrix.h>

/* constants for rendering tables */
const int PADDING = 1;
const int PRECISION = 6;
const char *FORMATTING  = "%.6f";

void setMatrixValues(double value, char type, Matrix matrix)
{
  for (int i=0;i<matrix->n;i++)
  {
    for (int j=0;j<matrix->m;j++)
    {
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
      case 'R':
        matrix->values[i][j] = (double)rand()/(double)(RAND_MAX/value);
        break;
      }
    }
  }
}

void copyMatrix(Matrix source, Matrix target)
{
  for (int i=0;i<source->n;i++)
  {
    for (int j=0;j<source->m;j++)
    {
      target->values[i][j] = source->values[i][j];
    }
  }
}

void copyColumn(Matrix source, int idx_s, Matrix target, int idx_t)
{
  for (int i=0;i<target->n;i++)
  {
    target->values[i][idx_t] = source->values[i][idx_s];
  }
}

void copyRow(Matrix source, int idx_s, Matrix target, int idx_t)
{
  for (int j=0;j<source->m;j++)
  {
    target->values[idx_t][j] = source->values[idx_s][j];
  }
}

void transposeMatrix(Matrix source, Matrix target)
{
  for (int i=0;i<target->n;i++)
  {
    for (int j=0;j<target->m;j++)
    {
      target->values[i][j] = source->values[j][i];
    }
  }
}

double maxValue(int m, double *array)
{
  int max = array[0];

  for (int i=1; i < m; i++)
  {
    if (max < array[i])
      max = array[i];
  }

  return max;
}

double matrixMax(Matrix matrix)
{
  double max = maxValue(matrix->m, matrix->values[0]);

  for (int i=1; i < matrix->n; i++)
  {
    double curr = maxValue(matrix->m, matrix->values[i]);
    if (max < curr)
      max = curr;
  }

  return max;
}

double sumMatrix(Matrix matrix)
{
  double sum = 0;
  for (int i=0;i<matrix->n;i++)
    {
      for (int j=0;j<matrix->m;j++)
        {
          sum = matrix->values[i][j];
        }
    }
  return sum;
}

double meanMatrix(Matrix matrix)
{
  double sum = sumMatrix(matrix);
  return sum / (matrix->n * matrix->m);
}


int numDigits(int n)
{
  int add = 0;
  if (n < 0)
    n = n * -1;

  if (n < 1)
    n = 1;
  else
    n = floor(log10(n)) + 1;

  return n;
}

int numDigitsDouble(double x)
{
  int n = floor(x);

  n = numDigits(n);
  if (x < 0) n++;

  return n;
}

void drawMatrix(Matrix matrix)
{
  double matrix_max = matrixMax(matrix);
  int max_width = numDigitsDouble(matrix_max);
  max_width = max_width + PADDING + PRECISION - 1;

   /* header */
  printf("\n");
  printf("(%d, %d) ", matrix->n, matrix->m);
  int offset = numDigits(matrix->n) + numDigits(matrix->m) + 5;

  for (int i=offset; i<(max_width*matrix->m); i++)
  {
    printf("-");
  }

  for (int i=0; i<matrix->n; i++)
  {
    printf("\n");
    for (int j=0; j<PADDING; j++)
    {
      printf("\n");
    }
    for (int j=0; j<matrix->m; j++)
    {
      double value = matrix->values[i][j];
      printf(FORMATTING, value);
      int blank_fill = max_width - numDigitsDouble(value);
      for (int i=0; i<blank_fill; i++)
      {
        printf(" ");
      }
    }
  }

  printf("\n");
  for (int i=0; i<(max_width*matrix->m); i++)
  {
    printf("-");
  }
  for (int i=0; i<3; i++) printf("\n");
}

/* scaling operations */
void scaleColumn(Matrix matrix, int idx, double scalar)
{
  for (int i=0;i<matrix->n;i++)
  {
    matrix->values[i][idx] = (matrix->values[i][idx]) * scalar;
  }
}

void scaleMatrix(Matrix matrix, double scalar)
{
  for (int i=0;i<matrix->m;i++)
  {
    scaleColumn(matrix, i, scalar);
  }
}

void absMatrix(Matrix matrix)
{
  for (int i=0; i<matrix->n; i++)
  {
    for (int j=0; j<matrix->m; j++)
    {
      matrix->values[i][j] = fabs(matrix->values[i][j]);
    }
  }
}

/* addition and subtraction */
void addColumn(Matrix target, int idx1, Matrix source, int idx2)
{
  for (int i=0; i<target->n; i++)
    {
      target->values[i][idx1] = target->values[i][idx1] + source->values[i][idx2];
    }
}

void addMatrix(Matrix target, Matrix source)
{
  for (int i=0; i<target->m; i++)
    {
      addColumn(target, i, source, i);
    }
}


void subtractColumn(Matrix target, int idx1, Matrix source, int idx2)
{
  for (int i=0; i<target->n; i++)
  {
    target->values[i][idx1] = target->values[i][idx1] - source->values[i][idx2];
  }
}

void subtractMatrix(Matrix target, Matrix source) {
  for (int i=0; i<target->m; i++)
  {
    subtractColumn(target, i, source, i);
  }
}

/* dot product */
double dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2)
{
  double x = 0.0;

  switch (orient)
  {
  case 'R':
    assert(matrix1->m == matrix2->m);
    for (int i=0; i<matrix1->m; i++)
    {
      x = x + (matrix1->values[idx1][i] * matrix2->values[idx2][i]);
    }
    break;
  case 'C':
    assert(matrix1->n == matrix2->n);
    for (int i=0; i<matrix1->n; i++)
    {
      x = x + (matrix1->values[i][idx1] * matrix2->values[i][idx2]);
    }
    break;
  }

  return x;
}

double dotProductV(Matrix matrix1, Matrix matrix2)
{
  double x = dotProduct('C', matrix1, 0, matrix2, 0);
  return x;
}

/* norm and normalization */
double norm(char orient, Matrix matrix, int idx)
{
  double x = sqrt(dotProduct(orient, matrix, idx, matrix, idx));
  return x;
}

double normV(Matrix matrix) {
  double x = norm('C', matrix, 0);
  return x;
}

void normalizeColumn(Matrix matrix, int idx) {
  double _norm = norm('C', matrix, idx);
  if (_norm != 0) /* if norm is zero, zero vector */
    scaleColumn(matrix, idx, 1/_norm);
}

/* projection */
void project(Matrix source1, int idx1, Matrix source2, int idx2,
             Matrix target, int idx_t)
{
  copyColumn(source2, idx2, target, idx_t);
  double st_dot = dotProduct('C', source1, idx1, source2, idx2);
  double norm_squared = dotProduct('C', source2, idx2, source2, idx2);

  if (norm_squared == 0)
    return;

  st_dot = st_dot / norm_squared;
  scaleColumn(target, idx_t, st_dot);
}

/*  outer product of self */
void outerMatrix(Matrix source, int idx_s, Matrix target)
{
  assert(source->n == target->n);
  assert(target->n == target->m);

  for (int i=0; i<target->n; i++)
  {
    for (int j=0; j<target->m; j++)
    {
      target->values[i][j] = source->values[i][idx_s] * source->values[j][idx_s];
    }
  }
}


/* matrix multiplication */
void multiplyMatrices(Matrix source1, Matrix source2, Matrix target)
{
  assert(source1->m == source2->n);
  assert(source1->n == target->n);
  assert(source2->m == target->m);
  for (int i=0; i<source1->n; i++)
  {
    for (int j=0; j<source2->m; j++)
    {
      double value = 0;
      for (int k=0; k<source1->m; k++)
      {
        value = value + (source1->values[i][k] * source2->values[k][j]);
      }
      target->values[i][j] = value;
    }
  }
}

