/**
  ******************************************************************************
  * @file    devTestGroupCalibration.cpp
  * @author  @victorcroisfelt
  * 		 @brandaogbs
  * @version V2.0
  * @date    May 13, 2019.
  * @brief   Testing the calibration methods.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
/* STM32F4 */
#include <stm32f4xx_hal.h>

/* Standard C */
#include <stdio.h>
#include <math.h>

/* CppUTest */
#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

/* Test */
#include <test.h>

/* User's */
#include <matrix.h>
#include <magdata.h>
#include <twostep.h>
#include <mve.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/

/* Private function ----------------------------------------------------------*/

/* Number of samples (N) */
short N;

/* Scale factor (sf) */
float sf = 350.0f;

/* Group Tests ---------------------------------------------------------------*/
TEST_GROUP(SelfCalibTests)
{ 
	void setup()
	{
		setupDWT();
	} 
	
	void teardown()
	{
	
	}
}
; 

/* Tests ---------------------------------------------------------------------*/
TEST(SelfCalibTests,twoStep)
{
	printf("Two-step algorithm:\n");

	for (N = 100; N < 351; N = N+25)
	{
		/* Declaring raw data vectors */
		matrix *x,*y,*z;	

		/* Loading some data */
		setupRawData(&x,&y,&z,N);

		/* Preparing to save the results of the calibration procedure */
		matrix *xCal,*yCal,*zCal,*ParC;

		/* #1 - Calibrating*/
		if (N == 300)
		{
			printf("**Execution time for each stage (samples %i)**\n",N);
			twoStep(&xCal,&yCal,&zCal,&ParC,x,y,z,sf,1);	
		}

		/* #2 - Calibrating */
		printf("Execution time of the whole process (samples %i):\n",N);
		DWT->CYCCNT = 0; /* resetting the register used for time counting */

		/* Performing the two-step calibration process */
		twoStep(&xCal,&yCal,&zCal,&ParC,x,y,z,sf,0);

		/* Validating the time */
		timeValidator(DWT->CYCCNT,600000000);

		/* Freeing */
		freeMatrix(x);
		freeMatrix(y);
		freeMatrix(z);
		freeMatrix(xCal);
		freeMatrix(yCal);
		freeMatrix(zCal);
		freeMatrix(ParC);
	}
}

TEST(SelfCalibTests,MVE)
{
	printf("MVE algorithm:\n");
	
	for (N = 100; N < 351; N = N + 25)
	{
		/* Declaring raw data vectors */
		matrix *x,*y,*z;	

		/* Loading some data */
		setupRawData(&x,&y,&z,N);

		/* Preparing to save the results of the calibration procedure */
		matrix *xCal,*yCal,*zCal,*ParC,*ParP;

		/* Taking some intial parameters */
		dataTreatment(&ParP,x,y,z);

		/* #1 - Calibrating */
		if (N == 300)
		{
			printf("**Execution time for each stage (samples %i)**\n",N);
			mve(&xCal,&yCal,&zCal,&ParC,x,y,z,ParP,sf,1);
		}

		/* #2 - Calibrating */
		printf("Execution time of the whole process (samples %i):\n",N);

		DWT->CYCCNT = 0; /* resetting the register used for time counting */

		/* Performing the MVE calibration process */
		mve(&xCal,&yCal,&zCal,&ParC,x,y,z,ParP,sf,0);

		timeValidator(DWT->CYCCNT,600000000);

		/* Freeing */
		freeMatrix(x);
		freeMatrix(y);
		freeMatrix(z);
		freeMatrix(ParP);
		freeMatrix(xCal);
		freeMatrix(yCal);
		freeMatrix(zCal);
		freeMatrix(ParC);
	}
}
/******************************* END OF FILE **********************************/