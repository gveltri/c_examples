/*
  @file matrix.h
  @author Gerardo Veltri
  Matrix manipulation
*/
#ifndef MATRIX_HEADER
#define MATRIX_HEADER

void setMatrixValues(float value, char type, Matrix matrix);
void copyMatrix(Matrix source, Matrix target);
void transposeMatrix(Matrix source, Matrix target);

void scaleColumn(Matrix matrix, int idx, float scalar);
void scaleMatrix(Matrix matrix, float scalar);
void absMatrix(Matrix matrix);

void addColumn(Matrix target, int idx1, Matrix source, int idx2);
void addMatrix(Matrix target, Matrix source);
void subtractColumn(Matrix target, int idx1, Matrix Source, int idx2);
void subtractMatrix(Matrix target, Matrix source);

void draw2DMatrix(Matrix matrix);

float sumMatrix(Matrix matrix);
float meanMatrix(Matrix matrix);
float matrixMax(Matrix matrix);

float dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2);
float dotProductV(Matrix matrix1, Matrix matrix2);

float norm(char orient, Matrix matrix, int idx);
float normV(Matrix matrix);
void normalizeColumn(Matrix matrix, int idx);

void outerMatrix(Matrix source, int idx_s, Matrix target);

void multiplyMatrices(Matrix source1, Matrix source2, Matrix target);

void project(Matrix source1, int idx1, Matrix source2, int idx2,
             Matrix target, int idx_t);


#endif
