/**
  ******************************************************************************
  * @file	 mve.cpp
  * @author  @victorcroisfelt    
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   MVE functions.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
#include <hdr\mve.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/

/**
  ==============================================================================
  @add:			   ----- MVE Private Functions -----
  ==============================================================================
  */

void sharedComp(matrix **xsubx0,matrix **ysuby0,matrix **zsubz0,matrix *x,matrix *y,matrix *z,matrix *ParP);

void b1tVector(matrix **b1t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b2tVector(matrix **b2t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b3tVector(matrix **b3t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b4tVector(matrix **b4t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b5tVector(matrix **b5t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b6tVector(matrix **b6t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b7tVector(matrix **b7t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b8tVector(matrix **b8t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void b9tVector(matrix **b9t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);

void BMatrix(matrix **B,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);

void MVEvar(matrix **r,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP);
void residue(matrix **e,matrix *r,matrix *x,float sf);
void qVector(matrix **q,matrix *B,matrix *e,matrix *x);

void cPar(matrix **ParC,matrix *q,matrix *ParP);
void calVec1(matrix **xCal,matrix **yCal,matrix **zCal,matrix *ParC,matrix *x,matrix *y,matrix *z);

/* Private functions ---------------------------------------------------------*/
/**
  ==============================================================================
  @def:			   ##### MVE Private Functions #####
  ==============================================================================
  */
/**
  * @brief Evalute the shared computations.
  */
void sharedComp(matrix **xsubx0,matrix **ysuby0,matrix **zsubz0,matrix *x,matrix * y,matrix *z,matrix *ParP)
{
	/* Colleting Offset Values */
	float *ptrParP = ParP->data;

	float x0 = ptrParP[0]; float y0 = ptrParP[1]; float z0 = ptrParP[2];

	/* 'x' - x0 */
	(*xsubx0) = makeMatrix(1,x->col);
	sumSubScalarMatrix(&(*xsubx0),x,-x0);	

	/* 'y' - y0 */
	(*ysuby0) = makeMatrix(1,y->col);
	sumSubScalarMatrix(&(*ysuby0),y,-y0);	

	/* 'z' - z0 */
	(*zsubz0) = makeMatrix(1,z->col);
	sumSubScalarMatrix(&(*zsubz0),z,-z0);
}

/**
  * @brief Evalute the b1 transposed vector.
  */
void b1tVector(matrix **b1t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

	/* b1 = (2 * sin(ro) * (a * (vy - yo) - b * sin(ro) * (vx - xo))) / (a^2 * b * cos(ro)^2) - (2 * vx - 2 * xo) / a^2 - (2 * (sin(la) * sin(ro) - ...
	 cos(la) * cos(ro) * sin(fi)) * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) * (vx - xo) + a * b * cos(ro) * (vz - zo) - a * c * sin(la) * (vy - yo))) ...
	 / (a^2 * b * c * cos(fi)^2 * cos(la)^2 * cos(ro)^2) */
	scaleMatrix(&aux0,xsubx0,-b*sin(ro));
	scaleMatrix(&aux1,ysuby0,a);
	sumMatrix(&aux2,aux1,aux0);
	scaleMatrix(&aux0,aux2,(2.0f*sin(ro))/(a*a*b*cos(ro)*cos(ro)));
	scaleMatrix(&aux1,xsubx0,-2.0f/(a*a));
	sumMatrix(&aux2,aux0,aux1);

	scaleMatrix(&aux0,zsubz0,a*b*cos(ro));
	scaleMatrix(&aux1,ysuby0,-a*c*sin(la));
	sumMatrix(&aux3,aux0,aux1);
	scaleMatrix(&aux0,xsubx0,b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi)));
	sumMatrix(&aux1,aux0,aux3);
	scaleMatrix(&aux0,aux1,-2.0f*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))/(a*a*b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	sumMatrix(&aux1,aux2,aux0);

	/* Attribution */
	transposeMatrix(&(*b1t),aux1);

	/* Freeing */
	freeMatrix(aux0); 
	freeMatrix(aux1); 
	freeMatrix(aux2);
	freeMatrix(aux3);
}

/**
  * @brief Evalute the b2 transposed vector.
  */
void b2tVector(matrix **b2t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

	/* b2 = (2 * sin(la) * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) * (vx - xo) + a * b * cos(ro) * (vz - zo) - a * c * sin(la) * (vy - yo))) / ...
	 (a * b^2 * c * cos(fi)^2 * cos(la)^2 * cos(ro)^2) - (2 * (a * (vy - yo) - b * sin(ro) * (vx - xo))) / (a * b^2 * cos(ro)^2) */
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	scaleMatrix(&aux1,zsubz0,a*b*cos(ro));
	sumMatrix(&aux2,aux1,aux0);
	scaleMatrix(&aux0,xsubx0,b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi)));
	sumMatrix(&aux1,aux0,aux2);
	scaleMatrix(&aux0,aux1,(2.0f*sin(la))/(a*b*b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	scaleMatrix(&aux1,ysuby0,a);
	scaleMatrix(&aux2,xsubx0,-b*sin(ro));
	sumMatrix(&aux3,aux1,aux2);
	scaleMatrix(&aux1,aux3,-2.0f/(a*b*b*cos(ro)*cos(ro)));

	sumMatrix(&aux2,aux0,aux1);

	/* Attribution */
	transposeMatrix(&(*b2t),aux2);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2);
	freeMatrix(aux3);
}

/**
  * @brief Evalute the b3 transposed vector.
  */
void b3tVector(matrix **b3t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);

	/* b3 = -(2 * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) * (vx - xo) + a * b * cos(ro) * (vz - zo) - a * c * sin(la) * (vy - yo))) ...
	 / (a * b * c^2 * cos(fi)^2 * cos(la)^2 * cos(ro)) */
	scaleMatrix(&aux0,xsubx0,(b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))));
	scaleMatrix(&aux1,zsubz0,a*b*cos(ro));
	sumMatrix(&aux2,aux0,aux1);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	sumMatrix(&aux1,aux2,aux0);

	scaleMatrix(&aux0,aux1,-2.0f/(a*b*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)));

	/* Attribution */
	transposeMatrix(&(*b3t),aux0);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
}

/**
  * @brief Evalute the b4 transposed vector.
  */
void b4tVector(matrix **b4t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);
	matrix *aux4 = makeMatrix(1,xsubx0->col);

	/* b4 = (2 * (a * (vy - yo) - b * sin(ro) * (vx - xo)) .* (vy - yo)) ./ (a^2 * b^2 * cos(ro)^2) ...
	 - (2 * (a * (vy - yo) - b * sin(ro) * (vx - xo)).^2) ./ (a^3 * b^2 * cos(ro)^2) ...
	 - (2 * (vx - xo).^2) ./ a^3 ...
	 - (2 * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) .* (vx - xo) + a * b * cos(ro) .* (vz - zo) - ...
	 a * c * sin(la) * (vy - yo)).^2) ./ (a^3 * b^2 * c^2 * cos(fi)^2 * cos(la)^2 * cos(ro)^2) ...
	 + (2 * (b * cos(ro) .* (vz - zo) - c * sin(la) * (vy - yo)) .* (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) .* (vx - xo) + a * b * cos(ro) .* (vz - zo) ... 
	 - a * c * sin(la) .* (vy - yo))) ./ (a^2 * b^2 * c^2 * cos(fi)^2 * cos(la)^2 * cos(ro)^2) */
	scaleMatrix(&aux0,ysuby0,a);
	scaleMatrix(&aux1,xsubx0,-b*sin(ro));
	sumMatrix(&aux2,aux0,aux1);
	mulElement(&aux0,aux2,ysuby0);
	scaleMatrix(&aux1,aux0,2.0f/(a*a*b*b*cos(ro)*cos(ro)));

	scaleMatrix(&aux0,ysuby0,a);
	scaleMatrix(&aux2,xsubx0,-b*sin(ro));
	sumMatrix(&aux3,aux0,aux2);
	squareElement(&aux0,aux3);
	scaleMatrix(&aux2,aux0,-2.0f/(a*a*a*b*b*cos(ro)*cos(ro)));

	sumMatrix(&aux3,aux1,aux2);

	squareElement(&aux0,xsubx0);
	scaleMatrix(&aux1,aux0,-2.0f/(a*a*a));

	sumMatrix(&aux2,aux3,aux1);

	scaleMatrix(&aux0,xsubx0,(b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))));
	scaleMatrix(&aux1,zsubz0,a*b*cos(ro));
	sumMatrix(&aux3,aux0,aux1);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	sumMatrix(&aux1,aux3,aux0);
	squareElement(&aux0,aux1);
	scaleMatrix(&aux1,aux0,-2.0f/(a*a*a*b*b*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	sumMatrix(&aux0,aux2,aux1);

	scaleMatrix(&aux1,xsubx0,(b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))));
	scaleMatrix(&aux2,zsubz0,a*b*cos(ro));
	sumMatrix(&aux3,aux1,aux2);
	scaleMatrix(&aux1,ysuby0,a*c*sin(la));
	sumMatrix(&aux2,aux3,aux1);
	scaleMatrix(&aux1,zsubz0,b*cos(ro));
	scaleMatrix(&aux3,ysuby0,-c*sin(la));
	sumMatrix(&aux4,aux1,aux3);
	mulElement(&aux1,aux4,aux2);
	scaleMatrix(&aux2,aux1,2.0f/(a*a*b*b*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	sumMatrix(&aux1,aux0,aux2);

	/* Attribution */
	transposeMatrix(&(*b4t),aux1);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
	freeMatrix(aux3);
	freeMatrix(aux4);
}

/**
  * @brief Evalute the b5 transposed vector.
  */
void b5tVector(matrix **b5t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);

	/* b5 =  -(2 * (vy - yo) .* (((a * c * sin(la)^2) * (vy - yo)) - ((a * b * cos(ro) * sin(la)) * (vz - zo)) - ((b * c * sin(la)^2 * sin(ro)) * (vx - xo)) + 
	 ((a * c * cos(fi)^2 * cos(la)^2) * (vy - yo)) - ((b * c * cos(fi)^2 * cos(la)^2 * sin(ro)) * (vx - xo)) + ((b * c * cos(la) * cos(ro) * sin(fi) * sin(la)) * (vx - xo))))
	 /(a * b^3 * c * cos(fi)^2 * cos(la)^2 * cos(ro)^2)  */
	scaleMatrix(&aux0,xsubx0,(a*c*sin(la)*sin(la)));
	scaleMatrix(&aux1,zsubz0,-(a*b*cos(ro)*sin(la)));
	sumMatrix(&aux2,aux0,aux1);
	scaleMatrix(&aux0,xsubx0,-(b*c*sin(la)*sin(la)*sin(ro)));
	sumMatrix(&aux1,aux2,aux0);
	scaleMatrix(&aux0,ysuby0,(a*c*cos(fi)*cos(fi)*cos(la)*cos(la)));
	sumMatrix(&aux2,aux1,aux0);
	scaleMatrix(&aux0,xsubx0,-(b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*sin(ro)));
	sumMatrix(&aux1,aux2,aux0);
	scaleMatrix(&aux0,xsubx0,(b*c*cos(la)*cos(ro)*sin(fi)*sin(la)));	
	sumMatrix(&aux2,aux1,aux0);
	mulElement(&aux0,ysuby0,aux2);

	scaleMatrix(&aux1,aux0,-2.0f/(a*b*b*b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	/* Attribution */
	transposeMatrix(&(*b5t),aux1);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
}

/**
  * @brief Evalute the b6 transposed vector.
  */
void b6tVector(matrix **b6t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);

	/* b6 = -(2 * (vz - zo) .* ( ( (a * b * cos(ro)) * (vz - zo) ) - ( (a * c * sin(la)) * (vy - yo) ) + ( (b * c * sin(la) * sin(ro)) * (vx - xo) ) - ...
     ( (b * c * cos(la) * cos(ro) * sin(fi)) * (vx - xo) ) ) ) ./ (a * b * c^3 * cos(fi)^2 * cos(la)^2 * cos(ro)) */
	scaleMatrix(&aux0,zsubz0,(a*b*cos(ro)));
	scaleMatrix(&aux1,ysuby0,-(a*c*sin(la)));
	sumMatrix(&aux2,aux0,aux1);
	scaleMatrix(&aux0,xsubx0,(b*c*sin(la)*sin(ro)));
	sumMatrix(&aux1,aux2,aux0);
	scaleMatrix(&aux0,xsubx0,-(b*c*cos(la)*cos(ro)*sin(fi)));
	sumMatrix(&aux2,aux0,aux1);
	mulElement(&aux0,zsubz0,aux2);

	scaleMatrix(&aux1,aux0,-2.0f/(a*b*c*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)));

	/* Attribution */
	transposeMatrix(&(*b6t),aux1);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
}

/**
  * @brief Evalute the b7 transposed vector.
  */
void b7tVector(matrix **b7t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

    /* b7 = -(2 * ( b * (vx - xo) - (a * sin(ro)) * (vy - yo) ) .* ( (a * c * sin(la)^2) * (vy - yo) - ...
     (a * b * cos(ro) * sin(la)) * (vz - zo) - (b * c * sin(la)^2 * sin(ro)) * (vx - xo) +...
     (a * c * cos(fi)^2 * cos(la)^2) * (vy - yo) - (b * c * cos(fi)^2 * cos(la)^2 * sin(ro)) * (vx - xo) +...
     (b * c * cos(la) * cos(ro) * sin(fi) * sin(la)) * (vx - xo) ) ) / (a^2 * b^2 * c * cos(fi)^2 * cos(la)^2 * cos(ro)^3) */
	scaleMatrix(&aux0,ysuby0,(a*c*sin(la)*sin(la)));
	scaleMatrix(&aux1,zsubz0,-(a*b*cos(ro)*sin(la)));
	sumMatrix(&aux2,aux0,aux1);
	scaleMatrix(&aux0,xsubx0,-(b*c*sin(la)*sin(la)*sin(ro)));
	sumMatrix(&aux1,aux2,aux0);
	scaleMatrix(&aux0,ysuby0,(a*c*cos(fi)*cos(fi)*cos(la)*cos(la)));
	sumMatrix(&aux2,aux1,aux0);
	scaleMatrix(&aux0,xsubx0,-(b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*sin(ro)));
	sumMatrix(&aux1,aux2,aux0);
	scaleMatrix(&aux0,xsubx0,(b*c*cos(la)*cos(ro)*sin(fi)*sin(la)));

	sumMatrix(&aux2,aux1,aux0);

	scaleMatrix(&aux0,xsubx0,b);
	scaleMatrix(&aux1,ysuby0,-(a*sin(ro)));

	sumMatrix(&aux3,aux0,aux1);

	mulElement(&aux0,aux3,aux2);

	scaleMatrix(&aux1,aux0,-2/(a*a*b*b*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)*cos(ro)));

	/* Attribution */
	transposeMatrix(&(*b7t),aux1);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
	freeMatrix(aux3);
}

/**
  * @brief Evalutes the b8 transposed vector.
  */
void b8tVector(matrix **b8t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

	/* b8 = (2 * sin(fi) * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) .* (vx - xo) +...
     a * b * cos(ro) .* (vz - zo) - a * c * sin(la) .* (vy - yo)).^2) ./ ...
     (a^2 * b^2 * c^2 * cos(fi)^3 * cos(la)^2 * cos(ro)^2) - (2 .* (vx - xo) .* ...
     (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) .* (vx - xo) + a * b * cos(ro) .* (vz - zo) - ...
     a * c * sin(la) .* (vy - yo))) ./ (a^2 * b * c * cos(fi) * cos(la) * cos(ro)) */
	scaleMatrix(&aux0,xsubx0,b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi)));
	scaleMatrix(&aux1,zsubz0,a*b*cos(ro));
	sumMatrix(&aux2,aux0,aux1);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	sumMatrix(&aux1,aux2,aux0);
	squareElement(&aux2,aux1);
	scaleMatrix(&aux1,aux2,(2*sin(fi))/(a*a*b*b*c*c*cos(fi)*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	scaleMatrix(&aux0,xsubx0,(b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))));
	scaleMatrix(&aux2,zsubz0,a*b*cos(ro));
	sumMatrix(&aux3,aux0,aux2);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	sumMatrix(&aux2,aux3,aux0);
	mulElement(&aux0,xsubx0,aux2);
	scaleMatrix(&aux2,aux0,-2/(a*a*b*c*cos(fi)*cos(la)*cos(ro)));

	sumMatrix(&aux0,aux1,aux2);

	/* Attribution */
	transposeMatrix(&(*b8t),aux0);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
	freeMatrix(aux3);
}

/**
  * @brief Evalute the b9 transposed vector.
  */
void b9tVector(matrix **b9t,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	float *ptrParP = ParP->data;

	float a = ptrParP[3];  float b = ptrParP[4];  float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

	/* b9 = (2 * (b * c * (cos(la) * sin(ro) + cos(ro) * sin(fi) * sin(la)).* (vx - xo) - ...
	 a * c * cos(la).* (vy - yo)) .* (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)) .* (vx - xo) + ...
	 a * b * cos(ro).* (vz - zo) - a * c * sin(la).* (vy - yo))) ./ (a ^ 2 * b ^ 2 * c ^ 2 * cos(fi) ^ 2 * cos(la) ^ 2 * cos(ro) ^ 2) + ...
	 (2 * sin(la) * (b * c * (sin(la) * sin(ro) - cos(la) * cos(ro) * sin(fi)).* (vx - xo) + a * b * cos(ro).* (vz - zo) - ...
	 a * c * sin(la).* (vy - yo)) .^ 2) ./ (a ^ 2 * b ^ 2 * c ^ 2 * cos(fi) ^ 2 * cos(la) ^ 3 * cos(ro) ^ 2); */
	scaleMatrix(&aux0,xsubx0,(b*c*(cos(la)*sin(ro)+cos(ro)*sin(fi)*sin(la))));
	scaleMatrix(&aux1,ysuby0,-a*c*cos(la));
	sumMatrix(&aux2,aux0,aux1);

	scaleMatrix(&aux0,xsubx0,(b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi))));
	scaleMatrix(&aux1,zsubz0,a*b*cos(ro));
	sumMatrix(&aux3,aux0,aux1);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	sumMatrix(&aux1,aux3,aux0);

	mulElement(&aux0,aux2,aux1);
	scaleMatrix(&aux1,aux0,2/(a*a*b*b*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	scaleMatrix(&aux0,xsubx0,b*c*(sin(la)*sin(ro)-cos(la)*cos(ro)*sin(fi)));
	scaleMatrix(&aux2,zsubz0,a*b*cos(ro));
	sumMatrix(&aux3,aux0,aux2);
	scaleMatrix(&aux0,ysuby0,-a*c*sin(la));
	squareElement(&aux2,aux0);
	scaleMatrix(&aux0,aux2,(2*sin(la))/(a*a*b*b*c*c*cos(fi)*cos(fi)*cos(la)*cos(la)*cos(la)*cos(ro)*cos(ro)));

	sumMatrix(&aux2,aux1,aux0);

	/* Attribution */
	transposeMatrix(&(*b9t),aux2);

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1);
	freeMatrix(aux2); 
	freeMatrix(aux3);
}

/**
  * @brief Build the B matrix with the b'x't vectors.
  * @param B: pass the address of a pointer to a matrix structure (RETURN).
  * @param xsubx0: pass the address of a matrix structure.
  * @param ysuby0: pass the address of a matrix structure.
  * @param zsubz0: pass the address of a matrix structure.
  * @param ParP: pass the address of a matrix structure.
  * @retval None.
  */
void BMatrix(matrix **B,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	/* Calculating bxt Vectors */
	matrix *b1t = makeMatrix(xsubx0->col,1);
	b1tVector(&b1t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b2t = makeMatrix(xsubx0->col,1);
	b2tVector(&b2t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b3t = makeMatrix(xsubx0->col,1);
	b3tVector(&b3t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b4t = makeMatrix(xsubx0->col,1);
	b4tVector(&b4t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b5t = makeMatrix(xsubx0->col,1);
	b5tVector(&b5t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b6t = makeMatrix(xsubx0->col,1);
	b6tVector(&b6t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b7t = makeMatrix(xsubx0->col,1);
	b7tVector(&b7t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b8t = makeMatrix(xsubx0->col,1);
	b8tVector(&b8t,xsubx0,ysuby0,zsubz0,ParP);

	matrix *b9t = makeMatrix(xsubx0->col,1);
	b9tVector(&b9t,xsubx0,ysuby0,zsubz0,ParP);

	/* Creating pointers to bxt and B datas */
	float *ptrb1t = b1t->data, *ptrb2t = b2t->data, *ptrb3t = b3t->data; 
	float *ptrb4t = b4t->data, *ptrb5t = b5t->data, *ptrb6t = b6t->data;
	float *ptrb7t = b7t->data, *ptrb8t = b8t->data, *ptrb9t = b9t->data;

	float *ptrB = (*B)->data;

	/* Building process */
	for (short i = 0; i < (*B)->row; i++)
	{
		for (short j = 0; j < (*B)->col; j++)
		{
			if (i == 0)
				ptrB[i * (*B)->col + j] = ptrb1t[j];
			
			if (i == 1)
				ptrB[i * (*B)->col + j] = ptrb2t[j];
			
			if (i == 2)
				ptrB[i * (*B)->col + j] = ptrb3t[j];
			
			if (i == 3)
				ptrB[i * (*B)->col + j] = ptrb4t[j];
			
			if (i == 4)
				ptrB[i * (*B)->col + j] = ptrb5t[j];
		
			if (i == 5)
				ptrB[i * (*B)->col + j] = ptrb6t[j];
			
			if (i == 6)
				ptrB[i * (*B)->col + j] = ptrb7t[j];
			
			if (i == 7)
				ptrB[i * (*B)->col + j] = ptrb8t[j];
			
			if (i == 8)
				ptrB[i * (*B)->col + j] = ptrb9t[j];
		}
	}	

	/* Freeing */
	freeMatrix(b1t);
	freeMatrix(b2t);
	freeMatrix(b3t);
	freeMatrix(b4t);
	freeMatrix(b5t);
	freeMatrix(b6t);
	freeMatrix(b7t);
	freeMatrix(b8t);
	freeMatrix(b9t);
}

/**
  * @brief Calcute variable 'r'.
  * @param r: pass the address of a pointer to a matrix structure (RETURN).
  * @param xsubx0: pass the address of a matrix structure.
  * @param ysuby0: pass the address of a matrix structure.
  * @param zsubz0: pass the address of a matrix structure.
  * @param ParP: pass the address of a matrix structure.
  * @retval None.
  */
void MVEvar(matrix **r,matrix *xsubx0,matrix *ysuby0,matrix *zsubz0,matrix *ParP)
{
	/* Collecting parameters */
	float *ptrParP = ParP->data;

	float a = ptrParP[3]; float b = ptrParP[4]; float c = ptrParP[5];
	float ro = ptrParP[6]; float fi = ptrParP[7]; float la = ptrParP[8];

	/* Creating auxiliar matrix */
	matrix *aux0 = makeMatrix(1,xsubx0->col);
	matrix *aux1 = makeMatrix(1,xsubx0->col);
	matrix *aux2 = makeMatrix(1,xsubx0->col);
	matrix *aux3 = makeMatrix(1,xsubx0->col);

    /* r = ((vx - xo).^2 / a^2 + (-b * sin(ro) .* (vx - xo) + a .* (vy - yo)).^2 / a^2 / b^2 / cos(ro)^2 + (b * c * (sin(ro) * sin(la) - cos(ro) * sin(fi) * cos(la))...
	 .* (vx - xo) - a * c * sin(la) .* (vy - yo) + a * b * cos(ro) .* (vz - zo)).^2 / a^2 / b^2 / c^2 / cos(ro)^2 / cos(fi)^2 / cos(la)^2) */
	squareElement(&aux0,xsubx0);
	scaleMatrix(&aux3,aux0,(1/(a*a)));	
	scaleMatrix(&aux0,xsubx0,-b*sin(ro));
	scaleMatrix(&aux1,ysuby0,a);
	sumMatrix(&aux2,aux0,aux1);
	squareElement(&aux0,aux2);
	scaleMatrix(&aux1,aux0,1/(a*a)/(b*b)/(cos(ro)*cos(ro)));	
	sumMatrix(&aux2,aux3,aux1);
	scaleMatrix(&aux0,xsubx0,b*c*(sin(ro)*sin(la)-cos(ro)*sin(fi)*cos(la)));
	scaleMatrix(&aux1,ysuby0,-a*c*sin(la));
	sumMatrix(&aux3,aux0,aux1);
	scaleMatrix(&aux0,zsubz0,a*b*cos(ro));
	sumMatrix(&aux1,aux3,aux0);
	squareElement(&aux0,aux1);
	scaleMatrix(&aux1,aux0,1/(a*a)/(b*b)/(c*c)/(cos(ro)*cos(ro))/(cos(fi)*cos(fi))/(cos(la)*cos(la)));
	sumMatrix(&aux0,aux2,aux1);

	/* Calculating 'rt' */
	transposeMatrix(&(*r),aux0);

	/* Freeing */
	freeMatrix(aux0); 
	freeMatrix(aux1); 
	freeMatrix(aux2);
	freeMatrix(aux3);
}

/**
  * @brief Computing the residue.
  */
void residue(matrix **e,matrix *r,matrix *x,float sf)
{
	/* Calculating 'e' (Nx1) */
	matrix *aux = makeMatrix(x->col,1);

	scaleMatrix(&aux,r,-1.0f);

	sumSubScalarMatrix(&(*e),aux,sf*sf);

	/* Freeing */
	freeMatrix(aux);
}

  /**
  * @brief Computing the q vector.
  */
void qVector(matrix **q,matrix *B,matrix *e,matrix *x)
{
	/* #1TransposeMatrix (Nx9) */
	matrix *Bt = makeMatrix(x->col,9);
	transposeMatrix(&Bt,B);

	/* #1MultiplyMatrix (9x9) = (9xN) * (Nx9) */
	matrix *mul0 = makeMatrix(9, 9);
	multiplyMatrix(&mul0, B, Bt);

	/* #1InverseMatrix (9x9) */
	matrix *inv = makeMatrix(9,9);
	inverseMatrix(&inv,mul0);

	/* #1MultiplyMatrix (9xN) = (9x9)*(9xN) */
	matrix *mul1 = makeMatrix(9,x->col);
	multiplyMatrix(&mul1,inv,B);

	/* #1MultiplyMatrix (9x1) = (9xN)*(Nx1) */
	multiplyMatrix(&(*q),mul1,e);

	/* Freeing */
	freeMatrix(Bt);
	freeMatrix(mul0);
	freeMatrix(inv);
	freeMatrix(mul1);
}

/**
  * @brief Calibration parameters.
  */
void cPar(matrix **ParC,matrix *q,matrix *ParP)
{
	/* Updating and Attribution of 'ParC' */
	sumMatrix(&(*ParC),ParP,q);
}

/**
  * @brief Comoute the calibrated vectors.
  * @param xCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param yCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param zCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param ParC: pass the address of a matrix structure.
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @retval None.
  */
void calVec1(matrix **xCal,matrix **yCal,matrix **zCal,matrix *ParC,matrix *x,matrix *y,matrix *z)
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
	scaleMatrix(&(*xCal),aux0,1.0f/pParC[3]);	

	/* 'yCal' =  (('y' - y0) / b - 'xCal' * sin(rho)) / cos(rho) */

	/* #1SumSubScalar: 'aux0' = 'yN' - y0 */
	sumSubScalarMatrix(&aux0,y,-pParC[1]);
		
	/* #1ScaleMatrix: 'aux1' = 'aux0' / b */
	scaleMatrix(&aux1,aux0,1.0f/pParC[4]);
		
	/* #2ScaleMatrix: 'aux0' = 'xCal' * sin(rho) */
	scaleMatrix(&aux0,*xCal,sin(pParC[6]));
		
	/* #1SubMatrix: 'aux2' = 'aux1' - 'aux0' */
	subMatrix(&aux2,aux1,aux0);	
		
	/* #3ScaleMatrix: 'yCal' = 'aux2' / cos(rho) */
	scaleMatrix(&(*yCal),aux2,1.0f/cos(pParC[6]));
	
	/* 'zCal' = (('z' - z0) / c - 'xCal' * sin(phi) * cos(lambda) - 'yCal' * sin(lambda) * cos(phi)) / cos(phi) * cos(lambda) */
	
	/* #1SumSubScalar: 'aux0' = 'zN' - z0 */
	sumSubScalarMatrix(&aux0,z,-pParC[2]);
	
	/* #1ScaleMatrix: 'aux1' = 'aux0' / c */
	scaleMatrix(&aux1,aux0,1.0f/pParC[5]);
	
	/* #2ScaleMatrix: 'aux0' = 'xN' * sin(phi) * cos(lambda) */
	scaleMatrix(&aux0,*xCal,sin(pParC[7])*cos(pParC[8]));
	
	/* #1SubMatrix: 'aux2' = 'aux1' - 'aux0' */
	subMatrix(&aux2,aux1,aux0);
	
	/* #3ScaleMatrix: 'aux0' = 'yN' * sin(lambda) * cos(phi) */
	scaleMatrix(&aux0,*yCal,sin(pParC[8])*cos(pParC[7]));
	
	/* #2SubMatrix: 'aux1' = 'aux2' - 'aux0' */
	subMatrix(&aux1,aux2,aux0);
	
	/* #4ScaleMatrix: 'zCal' = 'aux1' / (cos(phi) * cos(lambda)) */
	scaleMatrix(&(*zCal),aux1,1.0f/(cos(pParC[7])*cos(pParC[8])));	

	/* Freeing */
	freeMatrix(aux0);
	freeMatrix(aux1); 
	freeMatrix(aux2);
}

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @def:			 ##### MVE Exported Functions #####
  ==============================================================================
  */
/**
  * @brief Treats the raw data to evalutes the ParP (vector w/ the initial calibration parameters).
  * @param ParP: pass the address of a pointer to a matrix structure (RETURN).
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @retval None.
  */
void dataTreatment(matrix **ParP,
	matrix *x,
	matrix *y,
	matrix *z)
{
	/* Calculating max. and min. values of the RawData vectors */
	float x_max = 0.0f,x_min = 0.0f;

	/* #1MaxMin */
	maxMinVector(&x_max,&x_min,x); 
	
	float y_max = 0.0f,y_min = 0.0f;

	/* #2MaxMin */
	maxMinVector(&y_max,&y_min,y); 
	
	float z_max = 0.0f,z_min = 0.0f;

	/* #3MaxMin */
	maxMinVector(&z_max,&z_min,z);
	
	/* Initial offset */
	float x0 = (x_max+x_min)/2;
	float y0 = (y_max+y_min)/2;
	float z0 = (z_max+z_min)/2;
	
	/* Initial scale factor */
	float a = 1.0f,b = 1.0f,c = 1.0f;
	
	/* Initial misalignment angles */
	float rho = 0.0f,phi = 0.0f,lambda = 0.0f;
	
	/* Creating 'ParP' (9x1) */
	(*ParP) = makeMatrix(9,1);	

	/* Attribution of 'ParP' (9x1) */
	float *pParP = (*ParP)->data;
	pParP[0] = x0;  pParP[1] = y0;  pParP[2] = z0;
	pParP[3] = a;   pParP[4] = b;   pParP[5] = c;
    pParP[6] = rho; pParP[7] = phi; pParP[8] = lambda;
}

/**
  * @brief MVE self-calibration algorithm.
  * @param xCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param yCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param zCal: pass the address of a pointer to a matrix structure (RETURN).
  * @param ParC: pass the address of a pointer to a matrix structure (RETURN).
  * @param x: pass the address of a matrix structure.
  * @param y: pass the address of a matrix structure.
  * @param z: pass the address of a matrix structure.
  * @param ParP: pass the address of a matrix structure.
  * @param sf: give the scale factor of raw data.
  * @retval None.
  */
void mve(matrix **xCal,
	matrix **yCal,
	matrix **zCal, 
	matrix **ParC,
	matrix *x,
	matrix *y,
	matrix *z,
	matrix *ParP,
	float sf,
	short test)
{
	if (1 == test) /* already performing the time tests */ 
	{
		/*--------------------------------------------------------------------*/
		printf("Shared computations:\n");

		DWT->CYCCNT = 0; 
		/* Computing the shared computations */
		matrix *xsubx0,*ysuby0,*zsubz0;

		sharedComp(&xsubx0,&ysuby0,&zsubz0,x,y,z,ParP);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Matrix B:\n");

		DWT->CYCCNT = 0; 
		/* Computing the B matrix */
		matrix *B = makeMatrix(9,x->col);

		BMatrix(&B,xsubx0,ysuby0,zsubz0,ParP);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("MVE variable:\n");

		DWT->CYCCNT = 0; 
		/* Computing the MVE variable */
		matrix *r = makeMatrix(x->col,1);

		MVEvar(&r,xsubx0,ysuby0,zsubz0,ParP);

		/* Freeing */
		freeMatrix(xsubx0);
		freeMatrix(ysuby0);
		freeMatrix(zsubz0);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Residue:\n");

		DWT->CYCCNT = 0; 
		/* Computing the residue */
		matrix *e = makeMatrix(x->col,1);

		residue(&e,r,x,sf);

		/* Freeing */
		freeMatrix(r);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Vector 'q':\n");

		DWT->CYCCNT = 0; 
		/* Computing the qVector */
		matrix *q = makeMatrix(9,1);

		qVector(&q,B,e,x);

		/* Freeing */
		freeMatrix(B);
		freeMatrix(e);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Calibration parameters:\n");

		DWT->CYCCNT = 0; 
		/* Computing the calibration parameters */
		(*ParC) = makeMatrix(9,1);

		cPar(&(*ParC),q,ParP);

		timeValidator(DWT->CYCCNT,600000000);

		/*--------------------------------------------------------------------*/
		printf("Calibrated vectors:\n");

		DWT->CYCCNT = 0; 
		/* Computing the calibrated vectors */
		calVec1(&(*xCal),&(*yCal),&(*zCal),(*ParC),x,y,z);

		timeValidator(DWT->CYCCNT,600000000);
	}
	else
	{
		/* Computing the shared computations */
		matrix *xsubx0,*ysuby0,*zsubz0;

		sharedComp(&xsubx0,&ysuby0,&zsubz0,x,y,z,ParP);

		/*--------------------------------------------------------------------*/
		/* Computing the B matrix */
		matrix *B = makeMatrix(9,x->col);

		BMatrix(&B,xsubx0,ysuby0,zsubz0,ParP);

		/*--------------------------------------------------------------------*/
		/* Computing the MVEvariable */
		matrix *r = makeMatrix(x->col,1);

		MVEvar(&r,xsubx0,ysuby0,zsubz0,ParP);

		/* Freeing */
		freeMatrix(xsubx0);
		freeMatrix(ysuby0);
		freeMatrix(zsubz0);

		/*--------------------------------------------------------------------*/
		/* Computing the residue */
		matrix *e = makeMatrix(x->col,1);

		residue(&e,r,x,sf);

		/* Freeing */
		freeMatrix(r);

		/*--------------------------------------------------------------------*/
		/* Computing the qVector */
		matrix *q = makeMatrix(9,1);

		qVector(&q,B,e,x);

		/* Freeing */
		freeMatrix(B);
		freeMatrix(e);

		/*--------------------------------------------------------------------*/
		/* Computing the calibration parameters */
		(*ParC) = makeMatrix(9,1);

		cPar(&(*ParC),q,ParP);

		/*--------------------------------------------------------------------*/
		/* Computing the calibrated vectors */
		calVec1(&(*xCal),&(*yCal),&(*zCal),(*ParC),x,y,z);
	}
}
/******************************* END OF FILE **********************************/