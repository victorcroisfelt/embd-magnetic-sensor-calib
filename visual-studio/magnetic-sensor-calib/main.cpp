/**
  ******************************************************************************
  * @file    main.cpp
  * @author  @victorcroisfelt
  *			 @brandaogbs
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Main program body.
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/
#include <stm32f4xx_hal.h>
#include <CppUTest/CommandLineTestRunner.h>

/* Private typedef -----------------------------------------------------------*/

/* Private define ------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"
#endif

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/

/* Private functions ---------------------------------------------------------*/
	
/* Initialization */
void SysTick_Handler(void)
{
	HAL_IncTick();
	HAL_SYSTICK_IRQHandler();
}

/* Run all tests? Folder:"...\dev\dev\Tests" */
void runAllTests()
{
	const char *p = "";
	CommandLineTestRunner::RunAllTests(0,&p);
}

/**
  * @brief  Main program
  * @param  None
  * @retval None
  */
int main(void)
{	
	runAllTests();
	
	return 0;
}
