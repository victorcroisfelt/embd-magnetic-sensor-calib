/**
  ******************************************************************************
  * @file    matrix.h
  * @author  @victorcroisfelt
  *			     @brandaogbs        
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Header file for linear algebra functions.
  ******************************************************************************
  */
/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MATRIX_H
#define __MATRIX_H

/* Includes ------------------------------------------------------------------*/
/* Standard C */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Exported types ------------------------------------------------------------*/
/** 
  * @brief Matrix structure definition  
  */ 
typedef struct _matrix 
{
	short	row;		/* Specifies the number of rows of a matrix.
                           This parameter can be any value from 0 to 65536. */
	short	col;		/* Specifies the number of columns of a matrix.
                           This parameter can be any value from 0 to 65536. */
	float	*data;		/* Specifies a pointer to an array with the matrix data.
						*/
} matrix;

/* Exported constants --------------------------------------------------------*/

/* Exported macro ------------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/** 
  * @comment: It is defined that functions return objects of matrix** type,  
  * while their inputs are of matrix* types.
  */
/**
  ==============================================================================
  @add:  			 ##### Matrix Exported Functions #####
  ==============================================================================
  */
/**
  ==============================================================================
  @add:				    ---- Tools and Creation ----
  ==============================================================================
  */
matrix* makeMatrix(short row,short col); /* Create a matrix object. */
void    freeMatrix(matrix *m); /* Free the memory allocated for a matrix. */
void    printMatrix(matrix *m);	/* Print the matrix. */
void	  struct2matrix(float *p,short *rowp,short *colp,matrix *a); /* Convert a matrix to a double array in C. */
matrix* matrix2struct(float *A,short rowp,short colp); /* Convert a double array in C to a matrix. */
void  	elementMatrix(float *res,matrix *m,short r_el,short col_el); /* Take a desired element from a matrix. */
void	  isAllPositiveElements(bool *test,matrix *m); /* Verify if the elements of a matrix are all positive. */
void	  maxMinVector(float *max,float *min,matrix *a); /* Find the max/min values of a vector. */
void	  oneMatrix(matrix **one); /* Create a matrix of 1s. */
void	  eyeMatrix(matrix **eye); /* Create an identity (square) matrix. */

/**
  ==============================================================================
  @add:    			        ---- Operations ----
  ==============================================================================
  */
void	transposeMatrix(matrix **mt,matrix *m); /* Transpose a matrix. */
void	detMatrix(float *det,matrix *m); /* Compute the determinant of a matrix. */
void	traceMatrix(float *trace,matrix *m); /* Compute the trace of a matrix. */
void	mulElement(matrix **mulel,matrix *a,matrix *b); /* Multiplication among the elements of a matrix. */
void	inverseElement(matrix **inv,matrix *m); /* Compute the inverse of each element within a matrix. */
void	squareElement(matrix **sqr,matrix *m); /* Square the each element of a matrix. */
void	sumSubScalarMatrix(matrix **addSub,matrix *m,float scalar); /* Scalar sum. */
void	scaleMatrix(matrix **scale,matrix *m,float scalar); /* Scalar multiplication. */
void	sumMatrix(matrix **sum,matrix *a,matrix *b); /* Add matrix A and B (same size). */
void	subMatrix(matrix **sub,matrix *a,matrix *b); /* Subtract matrix A and B (same size). */
void	multiplyMatrix(matrix **mul,matrix *a,matrix *b); /* Multiply matrix A by matrix B (multiplication rule). */
void	inverseMatrix(matrix **inv,matrix *a); /* Compute the inverse of a matrix. */
void	dotMatrix(matrix **dot,matrix *a,matrix *b); /* Dot product between matrix A and B. */
void	crossMatrix3x3(matrix **cross,matrix *a,matrix *b); /* Cross porduct between matrix A and B both with a 3x3 size. */
void	crossVector(matrix **cross,matrix *a,matrix *b); /* Cross product between vectors A and B with a 3x1 lenght. */
void	norm(float *norm,matrix *vec); /* Compute the norm of a vector. */
void	normVector(matrix **normalized,matrix *vec); /* Normalize a vector. */

/* Private types -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macros ------------------------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

#endif /* __MATRIX_H */
/******************************* END OF FILE **********************************/


