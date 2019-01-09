/*
  @file matrix.h
  @author Gerardo Veltri
  Matrix manipulation
*/
#ifndef MATRIX_HEADER
#define MATRIX_HEADER

double maccess(Matrix matrix, int i, int j);
void mset(Matrix matrix, int i, int j, double value);


void fillMatrix(double values[], Matrix matrix);
void setMatrixValues(double value, char type, Matrix matrix);
void copyMatrix(Matrix source, Matrix target);
void copyRow(Matrix source, int idx_s, Matrix target, int idx_t);
void switchRow(Matrix matrix, int row1, int row2);
void transposeMatrix(Matrix source, Matrix target);

void scaleColumn(Matrix matrix, int idx, double scalar);
void scaleRow(Matrix matrix, int idx, double scalar);
void scaleMatrix(Matrix matrix, double scalar);
void absMatrix(Matrix matrix);

void addColumn(Matrix target, int idx1, Matrix source, int idx2);
void addMatrix(Matrix target, Matrix source);
void subtractColumn(Matrix target, int idx1, Matrix Source, int idx2);
void subtractRow(Matrix target, int idx1, Matrix source, int idx2);
void subtractMatrix(Matrix target, Matrix source);

void addRowScalarMultiple(Matrix target, int idx_t, double scalar, Matrix source, int idx_s);

void drawMatrix(Matrix matrix);

double sumMatrix(Matrix matrix, int _abs);
double meanMatrix(Matrix matrix, int _abs);
double matrixMax(Matrix matrix, int _abs);

double dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2);
double dotProductV(Matrix matrix1, Matrix matrix2);

double norm(char orient, Matrix matrix, int idx);
double normV(Matrix matrix);
void normalizeColumn(Matrix matrix, int idx);

void simpleProject(Matrix source1, int idx1, Matrix source2, int idx2,
		   Matrix target, int idx_t);
void project(Matrix source1, int idx1, Matrix source2, int idx2, double proj_scalar,
             Matrix target, int idx_t, double tscalar);

void outerMatrix(Matrix source, int idx_s, Matrix target);

void simpleMultiplyMatrices(Matrix source1, Matrix source2, Matrix target);
void multiplyMatrices(Matrix source1, int transpose1, Matrix source2, int transpose2,
		      Matrix target, double tscalar);

#endif
