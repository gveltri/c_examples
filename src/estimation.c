/*
  @file estimation.c
  @author Gerardo Veltri
  Estimations and approximations
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mem.h>
#include <matrix.h>
#include <factorization.h>

/*
  Ordinary Least Squares

  Using QR decomposition for speed and
  numerical stability

  Given:
  ------
  At * A * x = At * b
  A = QR

  (QR)t * QR * x = (QR)t * b
  Rt * Qt * Q * R * x = Rt * Qt * b
  Rt * R * x = Rt * Qt * b
  R * x = Qt * b
 
  @param A matrix of observations
  @param x coefficients of approximation
  @param b vector of values to be approximated
*/
void ordinaryLeastSquares(Matrix A, Matrix x, Matrix b)
{
	assert(A->n == b->n);
	assert(A->m == x->n);

	Matrix Qtb = allocMatrix(A->n, b->m);

	Matrix QR[] = {
		allocMatrix(A->n, A->n),
		allocMatrix(A->n, A->m),
	};

	gramSchmidtQR(A, QR, 0);

	multiplyMatrices(QR[0], 1, b, 0, Qtb, 0);

	drawMatrix(QR[1]);
	drawMatrix(x);
	drawMatrix(Qtb);

	backSubstitution(QR[1], x, Qtb);

	freeMatrix(QR[0]);
	freeMatrix(QR[1]);
	freeMatrix(Qtb);
}

/*
  Linear Regression

  Appends a column of ones to observations A
  and then executes ordinary least squares

  @param A matrix of observations
  @param x coefficients of approximation
  @param b vector of values to be approximated
*/
void linearRegression(Matrix A, Matrix x, Matrix b)
{
	Matrix _A = allocMatrix(A->n, A->m+1);
	copyMatrix(A, _A);
	for (int i=0; i<_A->n; i++)
	{
          mset(_A, i, _A->m-1, 1);
	}

	ordinaryLeastSquares(_A, x, b);
	freeMatrix(_A);
}
