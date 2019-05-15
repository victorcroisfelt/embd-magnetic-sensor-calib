# Numerical Error and Temporal Analysis of Embedded Magnetic Sensor Calibration using ARM Cortex-M4 Processor


## Abstract
This work describes the implementation and results related to the numerical error and the temporal analysis of the calibration methodology of an embedded magnetometer. The system tested is part of an attitude determination system for aerospace applications, which uses low-cost tri-axial MEMS (micro electromechanical systems) magnetometer and an ARM Cortex-M4 as the processing unit. In order to validate and develop an efficient methodology for in run calibration, two different calibration algorithms were evaluated: a Two-Step process and a Minimum Variance Estimator (MVE), both described in the literature. The sensor's model consists of nine parameters to be estimated: bias, scale factor and misalignment angles between axes. A Test Unit was developed to evaluate both the numerical and temporal results of calibration. It uses IEEE 754 floating point standard variables with 32-bit and 64-bit wide in the ARM Cortex-M4. The former can be executed within the processor's Float Point Unit and the later cannot. The numerical data are compared to the one obtained using a Matlab code simulation. Result shows that the Two-Step present numerical instability for 32-bit wide variables, while the Minimum Variance Estimator has an inferior temporal execution but did not present such numerical issues.

## Prerequisites

#### Software
* "[Microsoft Visual Studio Enterprise 2015 Version 14.0.25431.01 Update 3](https://docs.microsoft.com/en-us/visualstudio/releasenotes/vs2015-version-history)"
* "[Microsoft .NET Framework Version 4.7.03190](https://docs.microsoft.com/en-us/visualstudio/releasenotes/vs2015-version-history)"
* "[VisualGDB v5.2](https://visualgdb.com/history/)"

##### Visual GDB Packages
* STM32 Devices v2018.12R2
* OpenOCD v20170609-0.10.0
* Embedded Profiler and Fast Semihosting v3.1
* CppUTest v3.8r1

#### Hardware
"[STM32F4DISCOVERY](https://www.st.com/en/evaluation-tools/stm32f4discovery.html)"

## Content

### Matlab

### Visual Studio

## Running the Tests


## Citing this Repository and License
