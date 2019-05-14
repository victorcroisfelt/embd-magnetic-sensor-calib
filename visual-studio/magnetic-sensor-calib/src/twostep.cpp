/**
  ******************************************************************************
  * @file	 twostep.cpp
  * @author  @victorcroisfelt    
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Two-step functions.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
#include <hdr\twostep.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/
/**
  ==============================================================================
  @add:			 ----- Two-Step Private Functions -----
  ==============================================================================
  */
void sfNorm(matrix **xsf,matrix **ysf,matrix **zsf,matrix *x,matrix *y,matrix *z,float sf);
void zOffsetting(matrix **zoff,float *z_off,matrix *zsf);
void measureMatrix(matrix **A_X,matrix *xsf,matrix *ysf,matrix *zoff);
void measureVector(matrix **x_m,matrix *A_X,short sample_number);
void cPar(matrix **ParC,matrix *x_m,float sf,float z_off);
void calVec(matrix **xCal,matrix **yCal,matrix **zCal,matrix *ParC,matrix *x,matrix *y,matrix *z);

/* Private functions ---------------------------------------------------------*/
/**
  ==============================================================================
  @def:			  ----- Two-Step Private Functions -----
  ==============================================================================
  */
/**
  * @brief sf Data normalization.
  * @param xsf: pass the address of a pointer to a matrix structure (RETURN).
  * @param ysf: pass the address of a pointer to a matrix structure (RETURN).
  * @param zsf: pass the address of a pointer to a matrix structure (RETURN).
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @retval None.
  */
void sfNorm(matrix **xsf,matrix **ysf,matrix **zsf,matrix *x,matrix *y,matrix *z,float sf)
{
	(*xsf) = makeMatrix(1,x->col);	
	(*ysf) = makeMatrix(1,x->col);
	(*zsf) = makeMatrix(1,x->col);

	scaleMatrix(&(*xsf),x,1/sf);
	scaleMatrix(&(*ysf),y,1/sf);
	scaleMatrix(&(*zsf),z,1/sf);
}

/**
  * @brief Z axis offsetting.
  * @param zoff: pass the address of a pointer to a matrix structure (RETURN).
  * @param z_off: pass the address of a float variable (RETURN).
  * @param zsf: pass the address of a matrix structure.
  * @retval None.
  */
void zOffsetting(matrix **zoff,float *z_off,matrix *zsf)
{
	/* Evaluating max and min values of z coordinates */
	float z_max = 0.0f,z_min = 0.0f;

	maxMinVector(&z_max,&z_min,zsf);

	/* Calculating offset */
	(*z_off) = 100.0f*(z_max-z_min);
	
	/* Creating the a new z vector with offset performed (1xN) */
	(*zoff) = makeMatrix(1,zsf->col); 

	/* Offset: avoid points near zero */
	sumSubScalarMatrix(&(*zoff),zsf,(*z_off));
}

/**
  * @brief Construct the measure matrix.
  * @param A_X: pass the address of a pointer to a matrix structure (RETURN).
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param zoff: pass the address of a matrix structure.
  * @retval None.
  */
void measureMatrix(matrix **A_X,matrix *xsf, matrix *ysf,matrix *zoff)
{
	/* Creating the measure matrix 'A_X' (Nx9) */
	(*A_X) = makeMatrix(xsf->col,9);

	/* Pointers to matrix data */
	float *ptrA_X = (*A_X)->data,*ptrX = xsf->data,*ptrY = ysf->data,*ptrZo = zoff->data;
	
	/* Building process */
	for (short i = 0; i < (*A_X)->row; i++)
	{
		for (short j = 0; j < (*A_X)->col; j++)
		{
			if (j == 0)
				ptrA_X[i*(*A_X)->col+j]=(ptrX[i]*ptrX[i])/(ptrZo[i]*ptrZo[i]);
			
			if(j == 1)
				ptrA_X[i*(*A_X)->col+j]=(ptrX[i]*ptrY[i])/(ptrZo[i]*ptrZo[i]);
			
			if(j == 2)
				ptrA_X[i*(*A_X)->col+j]=ptrX[i]/ptrZo[i];
			
			if(j == 3)
				ptrA_X[i*(*A_X)->col+j]=(ptrY[i]*ptrY[i])/(ptrZo[i]*ptrZo[i]);
			
			if(j == 4)
				ptrA_X[i*(*A_X)->col+j]=ptrY[i]/ptrZo[i];
		
			if(j == 5)
				ptrA_X[i*(*A_X)->col+j]=ptrX[i]/(ptrZo[i]*ptrZo[i]);
			
			if(j == 6)
				ptrA_X[i*(*A_X)->col+j]=ptrY[i]/(ptrZo[i]*ptrZo[i]);
			
			if(j == 7)
				ptrA_X[i*(*A_X)->col+j]=1/ptrZo[i];
			
			if(j == 8)
				ptrA_X[i*(*A_X)->col+j]=1/(ptrZo[i]*ptrZo[i]);
		}
	}	
}

/**
  * @brief Evaluates the measure vector.
  * @param x_m: pass the address of a pointer to a matrix structure (RETURN).
  * @param A_X: pass the address of a matrix structure.
  * @param sample_number: give the number of samples of the collected data.
  * @retval None.
  */
void measureVector(matrix **x_m,matrix *A_X,short sample_number)
{
	/* #1TranposeMatrix: 'AX_T' (9xN) */
	matrix *A_XT = makeMatrix(9,sample_number);
	transposeMatrix(&A_XT, A_X);
	
	/* #1MultiplyMatrix: 'mul0' (9x9) = 'A_XT' * 'A_X' */
	matrix *mul0 = makeMatrix(9,9);
	multiplyMatrix(&mul0,A_XT,A_X);
	
	/* #1InverseMatrix: 'inv' (9x9) = 'mul0'^(-1) */
	matrix *inv = makeMatrix(9,9);
	inverseMatrix(&inv,mul0);
	
	/* #2MultiplyMatrix: 'mul1' (9xN) = 'inv' * 'AX_T' */
	mul0 = makeMatrix(9,sample_number);
	multiplyMatrix(&mul0,inv,A_XT);
	
	/* Freeing */
	freeMatrix(A_XT); 
	freeMatrix(inv);
	
	/* #1Creating: 'ones' matrix (Nx1) */
	matrix *ones = makeMatrix(sample_number,1);
	oneMatrix(&ones);
	
	/* #3MultiplyMatrix: 'x_m' (9x1) = 'mul1' * 'ones' */
	(*x_m) = makeMatrix(9,1);
	multiplyMatrix(&(*x_m),mul0,ones);	

	/* Freeing */
	freeMatrix(ones); 
	freeMatrix(mul0);	
}

/**
  * @brief Evaluate the calibration parameters obtained from the Two-Step algorithm.
  * @param ParC: pass the address of a pointer to a matrix structure (RETURN).
  * @param x_m: pass the address of a matrix structure.
  * @param z_off: give the z_off desired.
  * @retval None.
  */
void cPar(matrix **ParC,matrix *x_m,float sf,float z_off)
{
	/* Colleting the data from the measure vector */
	float *ptrXm = x_m->data;
	
	float alpha = ptrXm[0];	float beta = ptrXm[1]; 		float gamma = ptrXm[2];
	float delta = ptrXm[3]; float epsilon = ptrXm[4];	float chi = ptrXm[5];
	float eta = ptrXm[6];   float iota = ptrXm[7];    	float kappa = ptrXm[8];
	
	/* d: Auxiliar parameter */
	float d = (-2*beta*beta)-(2*beta*gamma*epsilon)+(2*delta*gamma*gamma)+(2*alpha*epsilon*epsilon)+(8*alpha*delta);
	
	/* Calculating the calibration parameters */
	float lambda = (0.5*epsilon)/sqrt(-delta);
	float x0 = -(epsilon*epsilon*chi-2*beta*eta+4*delta*chi-beta*epsilon*iota+2*gamma*delta*iota-gamma*epsilon*eta)/d;
	float y0 = -(gamma*gamma*eta+4*alpha*eta-2*beta*chi-beta*gamma*iota+2*alpha*epsilon*iota-gamma*epsilon*chi)/d;
	float z0 = -(beta*beta*iota-beta*gamma*eta-4*alpha*delta*iota+2*alpha*epsilon*eta-beta*epsilon*chi+2*gamma*delta*chi)/d;
	float c  = (1/sf)*sqrt(kappa-alpha*x0*x0-beta*x0*y0-gamma*x0*z0-delta*y0*y0-epsilon*y0*z0+z0*z0);
	float b  = c/sqrt(-delta);
	float a  = sqrt(2)*c*c*(1-lambda*lambda)*sqrt(-1/(2*alpha*c*c*(1-lambda*lambda)*(1-lambda*lambda)-(beta*b*lambda+gamma*c*lambda*lambda)*(beta*b*lambda+gamma*c)));
	float rho = a*(beta*b+gamma*c*lambda)/(2*c*c*(1-lambda*lambda));
	float phi = a*(beta*b*lambda+gamma*c)/(2*c*c*(1-lambda*lambda));

	/* Updating and atributting  the calibration parameters to a 9x1 vector 'ParC' */
	(*ParC) = makeMatrix(9,1);	
	float *pParC = (*ParC)->data;
	
	pParC[0] = x0*sf; pParC[1] = y0*sf; pParC[2] = (z0-z_off)*sf;
	pParC[3] = a*sf;  pParC[4] = b*sf;  pParC[5] = c*sf;
	pParC[6] = rho;	  pParC[7] = phi;   pParC[8] = lambda;
}

/**
  * @brief Calculate the calibrated vectors.
  * @param xCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param yCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param zCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param ParC: pass the address of a matrix structure.
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @retval None.
  */
void calVec(matrix **xCal,matrix **yCal,matrix **zCal,matrix *ParC,matrix *x,matrix *y,matrix *z)
{
	/* Collecting the calibration parameters from 'ParC' */
	float *pParC = ParC->data;
	
	/* Creating the calibrated vectors (1 x N) */
	(*xCal) = makeMatrix(1,x->col);
	(*yCal) = makeMatrix(1,x->col);
	(*zCal) = makeMatrix(1,x->col);
	
	/* Creating auxiliar matrix to evalute cal. vectors */
	matrix *aux0 = makeMatrix(1,x->col);
	matrix *aux1 = makeMatrix(1,x->col);
	matrix *aux2 = makeMatrix(1,x->col);	
	
	/* Calibrating data using the computed calibration parameters */
	
	/* 'xCal' = ('x' - x0) / a */

	/* #1SumSubScalar: 'aux0' = 'x' - x0 */
	sumSubScalarMatrix(&aux0,x,-pParC[0]);
		
	/* #1ScaleMatrix: 'xCal' = 'aux0' / a */
	scaleMatrix(&(*xCal),aux0,1/pParC[3]);	

	/* 'yCal' =  (('y' - y0) / b - 'xCal' * sin(rho)) / cos(rho) */

	/* #1SumSubScalar: 'aux0' = 'yN' - y0 */
	sumSubScalarMatrix(&aux0,y,-pParC[1]);
		
	/* #1ScaleMatrix: 'aux1' = 'aux0' / b */
	scaleMatrix(&aux1,aux0,1/pParC[4]);
		
	/* #2ScaleMatrix: 'aux0' = 'xCal' * sin(rho) */
	scaleMatrix(&aux0,*xCal,sin(pParC[6]));
		
	/* #1SubMatrix: 'aux2' = 'aux1' - 'aux0' */
	subMatrix(&aux2,aux1,aux0);	
		
	/* #3ScaleMatrix: 'yCal' = 'aux2' / cos(rho) */
	scaleMatrix(&(*yCal),aux2,1/cos(pParC[6]));
	
	/* 'zCal' = (('z' - z0) / c - 'xCal' * sin(phi) * cos(lambda) - 'yCal' * sin(lambda) * cos(phi)) / cos(phi) * cos(lambda) */
	
	/* #1SumSubScalar: 'aux0' = 'zN' - z0 */
	sumSubScalarMatrix(&aux0,z,-pParC[2]);
	
	/* #1ScaleMatrix: 'aux1' = 'aux0' / c */
	scaleMatrix(&aux1,aux0,1/pParC[5]);
	
	/* #2ScaleMatrix: 'aux0' = 'xN' * sin(phi) * cos(lambda) */
	scaleMatrix(&aux0,*xCal,sin(pParC[7])*cos(pParC[8]));
	
	/* #1SubMatrix: 'aux2' = 'aux1' - 'aux0' */
	subMatrix(&aux2,aux1,aux0);
	
	/* #3ScaleMatrix: 'aux0' = 'yN' * sin(lambda) * cos(phi) */
	scaleMatrix(&aux0,*yCal,sin(pParC[8])*cos(pParC[7]));
	
	/* #2SubMatrix: 'aux1' = 'aux2' - 'aux0' */
	subMatrix(&aux1,aux2,aux0);
	
	/* #4ScaleMatrix: 'zCal' = 'aux1' / (cos(phi) * cos(lambda)) */
	scaleMatrix(&(*zCal),aux1,1/(cos(pParC[7])*cos(pParC[8])));	

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1); 
	freeMatrix(aux2);
}
/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @def:			 ##### Two-Step Exported Functions #####
  ==============================================================================
  */
/**
  * @brief Two-Step self-calibration algorithm.
  * @param xCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param yCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param zCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param ParC: pass the address of a pointer to a matrix structure (RETURN).
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @param sf: give the scale factor of raw data.
  * @retval None.
  */
void twoStep(matrix **xCal,
	matrix **yCal,
	matrix **zCal, 
	matrix **ParC,
	matrix *x,
	matrix *y,
	matrix *z,
	float sf,
	short test)
{
	if (1 == test) /* already performing the time tests */  
	{
		/*--------------------------------------------------------------------*/
		printf("SF Normalization: \n");

		DWT->CYCCNT = 0; /* resetting the register used for time counting */
		/* Computing the sf normalization */
		matrix *xsf,*ysf,*zsf;

		sfNorm(&xsf,&ysf,&zsf,x,y,z,sf);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Z offset: \n");

		DWT->CYCCNT = 0; 
		/* Computing the offset over Z axis */
		matrix *zoff;
		float z_off = 0.0f;

		zOffsetting(&zoff,&z_off,zsf);

		freeMatrix(zsf);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Measure matrix: \n");

		DWT->CYCCNT = 0; /* resetting the register used for time counting */
		/* Computing of the measure matrix */
		matrix *A_X;

		measureMatrix(&A_X,xsf,ysf,zoff);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Measure vector: \n");

		DWT->CYCCNT = 0; /* resetting the register used for time counting */
		/* Computing the measure vector */
		matrix *x_m;

		measureVector(&x_m,A_X,x->col);

		freeMatrix(A_X);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Calibration parameters: \n");

		DWT->CYCCNT = 0; /* resetting the register used for time counting */
		/* Computing the calibration parameters */
		cPar(&(*ParC),x_m,sf,z_off);
	
		freeMatrix(x_m);
		freeMatrix(zoff);
		freeMatrix(xsf);
		freeMatrix(ysf);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Calibrated vectors: \n");

		DWT->CYCCNT = 0; /* resetting the register used for time counting */
		/* Computing the calibrated vectors */

		calVec(&(*xCal),&(*yCal),&(*zCal),(*ParC),x,y,z);

		timeValidator(DWT->CYCCNT,600000000);
	}
	else
	{
		/*--------------------------------------------------------------------*/
		/* Computing the sf normalization */
		matrix *xsf,*ysf,*zsf;

		sfNorm(&xsf,&ysf,&zsf,x,y,z,sf);
		/*--------------------------------------------------------------------*/
		/* Computing the offset over Z axis */
		matrix *zoff;
		float z_off = 0.0f;

		zOffsetting(&zoff,&z_off,zsf);

		freeMatrix(zsf);

		/*--------------------------------------------------------------------*/
		/* Computing of the measure matrix */
		matrix *A_X;

		measureMatrix(&A_X,xsf,ysf,zoff);

		/*--------------------------------------------------------------------*/
		/* Computing the measure vector */
		matrix *x_m;

		measureVector(&x_m,A_X,x->col);

		freeMatrix(A_X);

		/*--------------------------------------------------------------------*/
		/* Computing the calibration parameters */
		cPar(&(*ParC),x_m,sf,z_off);
	
		freeMatrix(x_m);
		freeMatrix(zoff);
		freeMatrix(xsf);
		freeMatrix(ysf);

		/*--------------------------------------------------------------------*/
		/* Computing the calibrated vectors */
		calVec(&(*xCal),&(*yCal),&(*zCal),(*ParC),x,y,z);
	}
}

/******************************* END OF FILE **********************************/