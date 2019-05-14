/**
  ******************************************************************************
  * @file    test.cpp
  * @author  @brandaogbs
  *          @victorcroisfelt
  * @version V2.0
  * @date    May 13, 2019.
  * @brief   Header files for functions to perform tests.
  ******************************************************************************
  */
/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __TEST_H
#define __TEST_H

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

/* CppUTest */
#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

/* Exported types ------------------------------------------------------------*/

/* Exported constants --------------------------------------------------------*/
#define DOUBLE_TOLERANCE 0.005

/* Exported macro ------------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @add:			##### Test Exported Functions #####
  ==============================================================================
  */
void setupDWT(void);
void timeValidator(unsigned cycles,unsigned timeout);
void checkMatrix(matrix *A,matrix *B);

/* Private types -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macros ------------------------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

#endif /* __TEST_H */
/******************************* END OF FILE **********************************/
