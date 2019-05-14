/**
  ******************************************************************************
  * @file    matrix.cpp
  * @author  @victorcroisfelt
  *			 @brandaogbs 
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Linear algebra functions.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
#include <hdr\matrix.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/
/**
  ==============================================================================
  @add:			      ##### Matrix Private Functions #####
  ==============================================================================
  */
void quicksort(matrix *a,short start,short end);
short partition(matrix *a,short start,short end);
void LUdecomposition(matrix *a,matrix **l,matrix **u);
matrix* solver(matrix *a,matrix *b);

/* Private functions ---------------------------------------------------------*/
/**
  ==============================================================================
  @def:			      ##### Matrix Private Functions #####
  ==============================================================================
  */
/**
  * @brief Perform the quicksort algorithm.
  * @param a: pass the address of a matrix struct (RETURN).
  * @param start: initial point (commonly 0).
  * @param end: end point (commonly vector's length).
  * @retval none.
  */
void quicksort(matrix *a,short start,short end) {
	short pivot;

	if (start < end) {
		pivot = partition(a,start,end);
		quicksort(a,start,pivot-1);
		quicksort(a,pivot+1,end);
	}
} 

/**
  * @brief This algorithm partitions an array so that all of the low values are
  * near the beginning and high values are near the end.
  * @param a: pass the address of a matrix struct (RETURN).
  * @param start: initial point.
  * @param end: end point.
  * @retval pivot: return the pivot encountered.
  */
short partition(matrix *a,short start,short end) {
	float *ptrA = a->data;
	float val = ptrA[start];
	short pivot = start;
	float temp;

	for (short i = start + 1; i <= end; i++) {
		if (ptrA[i] < val) {
			pivot++;
			temp = ptrA[pivot];
			ptrA[pivot] = ptrA[i];
			ptrA[i] = temp;
		}
	}

	temp = ptrA[start];
	ptrA[start] = ptrA[pivot];
	ptrA[pivot] = temp;

	return pivot;
}

/**
  * @brief Realize LU decomposition (A = L*U).
  * @param a: pass the address of a matrix struct.
  * @param l: pass the address of a pointer to a matrix struct (RETURN).
  * @param u: pass the address of a pointer to a matrix struct (RETURN).
  * @retval none.
  */
void LUdecomposition(matrix* a,matrix** l,matrix** u) {
	short i, j, k;
	float* ptrA;
	float* ptrL;
	float* ptrU;
	float sum;

	*l = makeMatrix(a->row,a->col);
	*u = makeMatrix(a->row,a->col);

	/* Step 1: Assign 1 to the diagonal of the lower matrix. */
	ptrL = (*l)->data;
	for (i = 0; i < a->col; i++) 
	{
		*ptrL = 1.0f;
		ptrL += a->col+1;
	}

	/* Step 2 */
	for (j = 0; j < a->col; j++) 
	{
	     /* Part A: Solve for the upper matrix. */
		for (i = 0; i <= j; i++) 
		{
			sum = 0.0f;
			for (k = 0; k < i; k++) 
				sum += (*l)->data[i*a->col+k]*(*u)->data[k*a->col+j];

			(*u)->data[i*a->col+j]=a->data[i*a->col+j]-sum;
		}
		/* Part B: Solve for the lower matrix. */
		for (i = j + 1; i < a->col; i++) 
		{
			sum = 0.0f;
			for (k = 0; k < j; k++) 
				sum += (*l)->data[i*a->col+k]*(*u)->data[k*a->col+j];
		
			(*l)->data[i*a->col+j]=1.0f/(*u)->data[j*a->col+j]*(a->data[i*a->col+j]-sum);
		}
	}

	return;
}

/**
  * @brief Resolve a linear system.
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval return the solvered matrix.
  */
matrix* solver(matrix* a,matrix* b) {
	short i, j, k;
	float sum;
	matrix *l = NULL;
	matrix *u = NULL;
	matrix *y;
	matrix *x;
	float* row;

	LUdecomposition(a,&l,&u);
	y = makeMatrix(a->row,1);
	x = makeMatrix(b->row,b->col);

	for (k = 0; k < b->col; k++) 
	{
	    /* Perform backward subsitituion with L */
	    /* L * y = B_k */
		for (i = 0; i < a->row; i++) 
		{
			row = l->data+i*a->col;
			sum = 0.0f;
			for (j = 0; j < i; j++) 
				sum += y->data[j]*(*row++);
		
			y->data[i] = (b->data[i*b->col+k]-sum)/(*row);
		}

		/* Perform backward subsitituion again with U */
		/* U * x = y */
		for (i = a->row - 1; i >= 0; i--) 
		{
			row = u->data+i*a->col+(a->col-1);
			sum = 0.0f;
			for (j = a->col - 1; j > i; j--) 
				sum += x->data[j*b->col+k]*(*row--);
			
			x->data[i*b->col+k]=(y->data[i]-sum)/(*row);
		}
	}

	/* Freeing */
	freeMatrix(l);
	freeMatrix(u);
	freeMatrix(y);

	return x;
}
/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @def:				   ##### Matrix Exported Functions #####
  ==============================================================================
  */
/**
  ==============================================================================
  @def:					  ---- Tools and Creation ----
  ==============================================================================
  */
/**
  * @brief Create a matrix.
  * @param row: number of rows.
  * @param col: number of cols.
  * @retval return the created matrix.
  */
matrix* makeMatrix(short row, short col) 
{
	matrix	*out;

	if (row == 0 || col == 0)
		return NULL;

	out = (matrix*) malloc(sizeof(matrix));
	
	if (out == NULL)
		return NULL;
    
	out->row = row;
	out->col = col;
	out->data = (float*) malloc(sizeof(float)*row*col);
    
	if (out->data == NULL)
		return NULL;
    
	memset(out->data,0.0f,row*col*sizeof(float));
	
	return out;
}

/**
  * @brief Free the memory allocated to store the matrix.
  * @param m: pass the address of a matrix.
  * @retval none.
  */
void freeMatrix(matrix *m) 
{
	if (m != NULL) 
	{
		if (m->data != NULL) 
		{
			free(m->data);
			m->data = NULL;
		}

		free(m);
		m = NULL;
	}
}

/**
  * @brief Print the values of the matrix.
  * @param m: pass the address of a matrix struct.
  * @retval none.
  */
void printMatrix(matrix *m) 
{
	float *ptr = m->data;      

	printf("%d %d\n",m->row,m->col);

	for (short i = 0; i < m->row; i++) 
	{
		for (short j = 0; j < m->col; j++)
			printf("%9.10f ",(*ptr++));
		printf("\n");
	}
}

/**
  * @brief Convert a matrix in a visible form in the debbuger.
  * @param p: pass the address of an array (RETURN).
  * @param rowp: pass the address of a short variable (RETURN).
  * @param colp: pass the address of a short variable (RETURN).
  * @param a: pass the address of a matrix struct.
  * @retval none.
  */
void struct2matrix(float *p,short *rowp,short *colp,matrix *a)
{
	float *ptr = a->data;
    
	*rowp = a->row;
	*colp = a->col;
	
	for (short i = 0; i < (a->row * a->col); i++)
		p[i] = *(ptr++);
}

/**
  * @brief Convert a matrix in C in a matrix type.
  * @param p: pass the address of an array.
  * @param rowp: number of rows.
  * @param colp: number of cols.
  * @retval return the matrix.
  */
matrix* matrix2struct(float *A,short rowp,short colp)
{
	matrix *m = makeMatrix(rowp,colp);
	float *ptr = m->data;

	m->row = rowp;
	m->col = colp;

	for (short i = 0; i < (rowp * colp); i++)
		*(ptr++) = A[i];
	
	return m;
}

/**
  * @brief Get an element of a matrix.
  * @param res: pass the address of a float variable (RETURN).
  * @param m: pass the address of a matrix struct.
  * @param r_el: number of the row.
  * @param col_el: number of the col.
  * @retval none.
  */
void elementMatrix(float *res,matrix *m,short r_el,short col_el)
{	
	float	*ptrM = m->data;
	
	for (short i = 0; i < m->row; i++)
	{
		for (short j = 0; j < m->col; j++)
		{
			if (i == (r_el - 1) && j == (col_el - 1))
				*res = *ptrM;
			
			ptrM++;
		}
	}
}

/**
  * @brief Search for negative elements in the matrix.
  * @param test: pass the address of a bool variable (RETURN).
  * @param m: pass the address of a matrix struct.
  * @retval none.
  */
void isAllPositiveElements(bool *test,matrix *m)
{
	float *pM = m->data;

	*test = true;
	
	for (short i = 0; i < (m->row*m->col); i++)
	{
		if (*(pM++) < 0.0f)*test = false;
	}
}


/**
  * @brief Search for the max. and min. values of a vector.
  * @param max: pass the address of a float variable (RETURN).
  * @param min: pass the address of a float variable (RETURN).
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void maxMinVector(float *max,float *min,matrix *a)
{
	float *ptrA = a->data;

	/* Assume first element as maximum and minimum */
	*max = ptrA[0];
	*min = ptrA[0];

	/* Find maximum and minimum */
	for (short i = 1; i < a->col; i++)
	{
	    /* If current element is greater than max */
		if (ptrA[i] > (*max))
			*max = ptrA[i];

		/* If current element is smaller than min */
		if (ptrA[i] < (*min))
			*min = ptrA[i];
	}
}

/**
  * @brief Generate a matrix with all elements equal 1.
  * @param one: pass the address of a pointer to a matrix struct (RETURN).
  * @retval none.
  */
void oneMatrix(matrix **one)
{
	float *ptrOne = (*one)->data;
	for (short i = 0; i < ((*one)->row * (*one)->col); i++)
		*(ptrOne++) = 1;
}

/**
  * @brief Create a identity (square) matrix.
  * @param eye: pass the address of a pointer to a matrix struct (RETURN).
  * @retval none.
 **/
void eyeMatrix(matrix **eye) 
{
	float *ptr = (*eye)->data;

	for (short i = 0; i < (*eye)->row; i++) 
	{
		*ptr = 1.0f;
		ptr += (*eye)->row+1;
	}
}

/**
  ==============================================================================
  @def:						 ---- Operations ----
  ==============================================================================
  */
/**
  * @brief Transpose a matrix (A -> AT).
  * @param mt: pass the address of a pointer to a matrix struct (RETURN).
  * @param m: pass the address of a matrix struct.
  * @retval none.
  */
void transposeMatrix(matrix **mt,matrix *m) 
{
	float *ptrM = m->data, *ptrMT;

	for (short i = 0; i < m->row; i++) 
	{
		ptrMT = &(*mt)->data[i];
		for (short j = 0; j < m->col; j++) 
		{
			*ptrMT = *(ptrM++);
			ptrMT += (*mt)->col;
		}
	}
}

/**
  * @brief Compute the determinant of a matrix using the LU decomposition.
  * @param det: pass the address of a float variable (RETURN).
  * @param m: pass the address of a matrix struct.
  * @retval none.
  */
void detMatrix(float *det,matrix *m) 
{
	*det = 1.0f;
    matrix* l = NULL;
    matrix* u = NULL;

    LUdecomposition(m,&l,&u);

    for (short i = 0; i < m->col; i++)
        *det *= u->data[i * m->col + i];

    freeMatrix(l);
    freeMatrix(u);
}

/**
  * @brief Compute the trace of a matrix.
  * @param trace: pass the address of a float variable (RETURN).
  * @param m: pass the address of a matrix struct.
  * @retval none.
  */
void traceMatrix(float *trace,matrix *m) 
{
	short	size;
	float	*ptr = m->data;
	*trace = 0.0f;

	if (m->row < m->col) 
		size = m->row;
	else 
		size = m->col;

	for (short i = 0; i < size; i++) 
	{
		(*trace) += *ptr;
		ptr += m->col + 1;
	}
}

/**
  * @brief Multiply element by element of a matrix.
  * @param mulel: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void mulElement(matrix **mulel,matrix *a,matrix *b)
{
	float *ptrA = a->data,*ptrB = b->data,*ptrMulel = (*mulel)->data;
	
	for (short i = 0; i < (a->row*a->col); i++)
		(*ptrMulel++) = (*ptrA++)*(*ptrB++);
}

/**
  * @brief Compute the inverse of each element of a matrix.
  * @param inv: pass the address of a pointer to a matrix struct (RETURN).
  * @param m: pass the address of a matrix struct.  
  * @retval none.
  */
void inverseElement(matrix **inv,matrix *m)
{
	float	*ptrM = m->data,*ptrInv = (*inv)->data;
	
	for (short i = 0; i < (m->row*m->col); i++)
	{
		if (*ptrM != 0.0)
			(*ptrInv++) = 1.0f/(*ptrM++);	
	}
}

/**
  * @brief Square each element of a matrix.
  * @param sqr: pass the address of a pointer to a matrix struct (RETURN).
  * @param m: pass the address of a matrix struct.  
  * @retval none.
  */
void squareElement(matrix **sqr,matrix *m)
{
	float *ptrM = m->data,*ptrSqr = (*sqr)->data;
	
	for (short i = 0; i < (m->row*m->col); i++)
	{
		(*ptrSqr++) = (*ptrM)*(*ptrM);	
		ptrM++;	
	}
		
}

/**
  * @brief Scalar sum.
  * @param addSub: pass the address of a pointer to a matrix struct (RETURN).
  * @param m: pass the address of a matrix struct.  
  * @param scalar: give the scalar desired to sum with all elements of a matrix.  
  * @retval none.
  */
void sumSubScalarMatrix(matrix **addSub,matrix *m,float scalar) 
{
	float *ptrM = m->data, *ptrAddSub = (*addSub)->data;

	for (short i = 0; i < (m->col*m->row); i++) 
		(*ptrAddSub++) = (*ptrM++)+scalar;
}

/**
  * @brief Scalar multiplication.
  * @param scale: pass the address of a pointer to a matrix struct (RETURN).
  * @param m: pass the address of a matrix struct.  
  * @param scalar: give the scalar desired to mul by each element of a matrix.
  * @retval none.
  */
void scaleMatrix(matrix **scale,matrix *m,float scalar) 
{
	float *ptrM = m->data,*ptrScale = (*scale)->data;

	for (short i = 0; i < (m->col*m->row); i++) 
		(*ptrScale++) = (*ptrM++)*scalar;
}

/**
  * @brief Sum two matrices (A + B).
  * @param sum: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void sumMatrix(matrix **sum,matrix *a,matrix *b)
{	
	float *ptrA = a->data,*ptrB = b->data,*ptrSum = (*sum)->data;

	for (short i = 0; i < (a->row*a->col); i++) 
		(*ptrSum++) = (*ptrA++) + (*ptrB++); 
}

/**
  * @brief Subtract two matrices (A - B).
  * @param sub: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void subMatrix(matrix **sub,matrix *a,matrix *b)
{
	float *ptrA = a->data,*ptrB = b->data,*ptrSub = (*sub)->data;

	for (short i = 0; i < (a->row*a->col); i++) 
		(*ptrSub++) = (*ptrA++)-(*ptrB++); 
}

/**
  * @brief Multiply two matrices (A*B = AB).
  * @param mul: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
 **/
void multiplyMatrix(matrix **mul,matrix *a,matrix *b) 
{
	float *ptrA,*ptrB,*ptrMul = (*mul)->data;
	
	for (short i = 0; i < a->row; i++) 
	{
		for (short j = 0; j < b->col; j++) 
		{
			ptrA = &a->data[i*a->col];
			ptrB = &b->data[j];
			*ptrMul = 0;
			for (short k = 0; k < a->col; k++) 
			{
				(*ptrMul) += (*ptrA)*(*ptrB);
				ptrA++;
				ptrB += b->col;
			}
			ptrMul++;
		}
	}
}

/**
  * @brief Compute the inverse of a matrix by using the LU decomposition.
  * @param inv: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @retval none.
  */
void inverseMatrix(matrix **inv,matrix *a) 
{	
	matrix *eye = makeMatrix(a->row,a->col);

	eyeMatrix(&eye);
	
	*inv = solver(a, eye);

	freeMatrix(eye);
}

/**
  * @brief Compute the dot product.
  * @param dot: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void dotMatrix(matrix **dot,matrix *a,matrix *b) 
{	
	float *ptrOut, *ptrA, *ptrB;
	
	if (a->col == 1 && a->row == 1)
		*dot = makeMatrix(1,1);
	else
		*dot = makeMatrix(1,a->col);
	
	ptrOut = (*dot)->data;		
	
	for (short i = 0; i < a->col; i++) 
	{
		ptrA = &a->data[i];
		ptrB = &b->data[i];
		for (short j = 0; j < a->row; j++) 
			(*ptrOut) = (*ptrOut) + *(ptrA + a->col * j)*(*(ptrB + a->col * j));
		
		if (a->row == 1 || a->col == 1)
			;
		else
			ptrOut++;
	}
}

/**
  * @brief Compute the cross product of 3x3 matrices.
  * @param cross: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void crossMatrix3x3(matrix **cross,matrix *a,matrix *b)
{	
	float *ptrA = a->data,*ptrB = b->data,*ptrCross = (*cross)->data;
	
	/* c(1,  :) = a(2,  :).*b(3,  :) - a(3,  :).*b(2,  :);  */
	*(ptrCross+0) = *(ptrA+a->col+0) * *(ptrB+2*a->col+0) - *(ptrA+2*a->col+0) * *(ptrB+a->col+0);
	*(ptrCross+1) = *(ptrA+a->col+1) * *(ptrB+2*a->col+1) - *(ptrA+2*a->col+1) * *(ptrB+a->col+1);
	*(ptrCross+2) = *(ptrA+a->col+2) * *(ptrB+2*a->col+2) - *(ptrA+2*a->col+2) * *(ptrB+a->col+2);
	
	/* c(2,  :) = a(3,  :).*b(1,  :) - a(1,  :).*b(3,  :);  */	 
	*(ptrCross+a->col+0) = *(ptrA+2*a->col+0) * *(ptrB+0) - *(ptrA+0) * *(ptrB+2*a->col+0);
	*(ptrCross+a->col+1) = *(ptrA+2*a->col+1) * *(ptrB+1) - *(ptrA+1) * *(ptrB+2*a->col+1);	
	*(ptrCross+a->col+2) = *(ptrA+2*a->col+2) * *(ptrB+2) - *(ptrA+2) * *(ptrB+2*a->col+2);
	
	/* c(3,  :) = a(1,  :).*b(2,  :) - a(2,  :).*b(1,  :);  */
	*(ptrCross+2*a->col+0) = *(ptrA+0) * *(ptrB+a->col+0) - *(ptrA+a->col+0) * *(ptrB+0);
	*(ptrCross+2*a->col+1) = *(ptrA+1) * *(ptrB+a->col+1) - *(ptrA+a->col+1) * *(ptrB+1);
	*(ptrCross+2*a->col+2) = *(ptrA+2) * *(ptrB+a->col+2) - *(ptrA+a->col+2) * *(ptrB+2);
}

/**
  * @brief Compute the cross product of 3x1 vectors.
  * @param cross: pass the address of a pointer to a matrix struct (RETURN).
  * @param a: pass the address of a matrix struct.
  * @param b: pass the address of a matrix struct.
  * @retval none.
  */
void crossVector(matrix **cross,matrix *a,matrix *b)
{	 
	float *ptrA = a->data,*ptrB = b->data,*ptrCross = (*cross)->data;
	
	*(ptrCross+0) = *(ptrA+1) * *(ptrB+2) - *(ptrA+2) * *(ptrB+1);
	*(ptrCross+1) = *(ptrA+2) * *(ptrB+0) - *(ptrA+0) * *(ptrB+2);
	*(ptrCross+2) = *(ptrA+0) * *(ptrB+1) - *(ptrA+1) * *(ptrB+0);
}

/**
  * @brief Compute the norm of a vector.
  * @param norm: pass the address of a float variable (RETURN).
  * @param vec: pass the address of a matrix struct.
  * @retval none.
  */
void norm(float *norm,matrix *vec)
{
	float *ptrVec = vec->data,sum = 0.0f;

	for (short i = 0; i < (vec->row * vec->col); i++)
		sum = sum+(*ptrVec*(*ptrVec++));
	
	 *norm = sqrt(sum);	
}

/**
  * @brief Normalize a vector.
  * @param normalized: pass the address of a pointer to a matrix struct (RETURN).
  * @param vec: pass the address of a matrix struct.
  * @retval none.
  */
void normVector(matrix **normalized,matrix *vec)
{		
	float value = 0.0f;
	
	norm(&value,vec);
	
	scaleMatrix(&(*normalized),vec,1.0f/value);
}
/******************************* END OF FILE **********************************/



