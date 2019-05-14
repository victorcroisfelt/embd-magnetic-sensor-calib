/**
  ******************************************************************************
  * @file    test.cpp
  * @author  @brandaogbs
  *          @victorcroisfelt
  * @version V2.0
  * @date    May 13, 2019.
  * @brief   Functions to perform tests.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
#include <hdr\test.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @def:			 ##### Test Exported Functions #####
  ==============================================================================
  */
/**
  * @brief DWT setup.
  * @retval none.
  */
void setupDWT(void)
{
	CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
	DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}

/**
  * @brief timeValidator test.
  * @param cycles: cycles number.
  * @param timeout: time limit.
  * @retval none.
  */
void timeValidator(unsigned cycles,unsigned timeout)
{
	if (cycles < timeout)
		printf("%d < %d cycles \n",cycles,timeout);
	else
	{
		char sz[64];
		sprintf(sz,"Violation: %d > %d cycles",cycles,timeout);
		FAIL(sz);
	}
}

/**
  * @brief Compare two matrix given an error tolerance.
  * @retval none.
 **/
void checkMatrix(matrix *A,matrix *B)
{
	float *pA = A->data;
	float *pB = B->data;

	for (short i = 0; i < (A->row*A->col); i++)
	{
		volatile float tolerance = (*pA > 0) ? (*pA) * DOUBLE_TOLERANCE : (*pA) * -DOUBLE_TOLERANCE;
		DOUBLES_EQUAL(*(pA++),*(pB++),tolerance)
	}
}
/******************************* END OF FILE **********************************/
