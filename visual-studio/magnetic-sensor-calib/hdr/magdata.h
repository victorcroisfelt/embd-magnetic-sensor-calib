/**
  ******************************************************************************
  * @file    magdata.h
  * @author  @victorcroisfelt   
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Getting some random data acquired from a magnetometer inserted into
  a Helmholtz coil. 
  * See more in:
  * [1] D. S. Batista, F. Granziera, M. C. Tosin and L. F. de Melo, “Three-Axial
  * Helmholtz Coil Design and Validation for Aerospace Applications,” in IEEE 
  * Transactions on Aerospace and Electronic Systems, vol. 54, no. 1, pp. 
  * 392-403, Feb. 2018.
  ******************************************************************************
  */
/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MAGDATA_H
#define __MAGDATA_H

/* Includes ------------------------------------------------------------------*/
/* User's */
#include <hdr\matrix.h>

/* Exported types ------------------------------------------------------------*/

/* Exported constants --------------------------------------------------------*/

/* Exported macro ------------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @add:				##### MagData Exported Functions #####
  ==============================================================================
  */
void setupRawData(matrix **x,			
	matrix **y, 
	matrix **z,
	short sample_number);

/* Private types -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macros ------------------------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

#endif /* __MAGDATA_H */
/******************************* END OF FILE **********************************/
