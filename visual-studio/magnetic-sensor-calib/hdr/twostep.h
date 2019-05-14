/**
  ******************************************************************************
  * @file    twostep.h
  * @author  @victorcroisfelt  
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Header file for the two-step calibration procedure.
  ******************************************************************************
  */
/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __TWOSTEP_H
#define __TWOSTEP_H

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
  @add:			##### Two-Step Process Exported Functions #####
  ==============================================================================
  */
void twoStep(matrix **xCal,			/* Find the self-calibration parameters and calibrate the vectors w/ the two-step calibration procedure. */									 
	matrix **yCal,
	matrix **zCal,
	matrix **ParC,
	matrix *x,
	matrix *y,
	matrix *z, 
	float sf,
	short test);
/* Private types -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macros ------------------------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

#endif /* __TWOSTEP_H */
/******************************* END OF FILE **********************************/
