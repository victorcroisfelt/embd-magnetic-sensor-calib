/**
  ******************************************************************************
  * @file    mve.h
  * @author  @victorcroisfelt  
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Header file for the MVE calibration procedure.
  ******************************************************************************
  */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MVE_H
#define __MVE_H

/* Includes ------------------------------------------------------------------*/
/* Standard C */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Matrix */
#include <hdr\matrix.h>

/* STM32F4 */
#include <stm32f4xx_hal.h>

/* Test */
#include <hdr\test.h>

/* Exported types ------------------------------------------------------------*/

/* Exported constants --------------------------------------------------------*/

/* Exported macro ------------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @add:			##### SelfCalibration Exported Functions #####
  ==============================================================================
  */
void dataTreatment(matrix **ParP, /* Treat the raw data to find the initial parameters (ParP) for calibration. */
	matrix *x,
	matrix *y,
	matrix *z);	
void mve(matrix **xCal, /* Find the self-calibration parameters and calibrate the vectors with the MVE calibration method. */										 
	matrix **yCal,
	matrix **zCal, 
	matrix **ParC,
	matrix *x,
	matrix *y,
	matrix *z,
	matrix *ParP,
	float sf,
	short test);
/* Private types -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macros ------------------------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

#endif /* __MVE_H */
/******************************* END OF FILE **********************************/




