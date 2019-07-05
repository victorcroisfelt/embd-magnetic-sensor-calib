# Numerical Error and Temporal Analysis of Embedded Magnetic Sensor Calibration using ARM Cortex-M4 Processor

This is a research-oriented code package that is primarily intended to allow readers to replicate the results of the article mentioned below and also encourage and accelerate further research on this topic:

Victor Croisfelt Rodrigues, Guilherme Brandão da Silva, Daniel Santana Batista, Marcelo Carvalho Tosin, and Francisco Granziera Jr, "[Numerical Error and Temporal Analysis of Embedded Magnetic Sensor Calibration using ARM Cortex-M4 Processor](https://doi.org/10.1109%2Finduscon.2018.8627205)," 13th IEEE INDUSCON, pp. 1387-1394, 2018.

The package is comprised of two parts: (a) some code written in C++/C that implements auto-calibration routines for embedded magnetic sensor data, and (b) a Matlab code that plot figures and simulates some data described in the temporal/numerical analysis part of the manuscript. To contextualize, in the sequel, we present the abstract of the article and other important information. Before proceeding, observe that the C++/C codes are available in a Visual Studio project which makes use of the VisualDGB package, allowing us to load the code via JTAG to a ARM Cortex-M4 microprocessor assembled in the STM32F4DISCOVERY developer kit. Plug-and-play support for other microprocessors is not provided here, but this is easy to get, since you only need to use the codes available in the header ("hrd") and source ("scr") folders (visualize them on "visual-studio/magnetic-sensor-calib/").

I hope this content helps in your reaseach and contributes to building the precepts behind open science. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions and if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## Abstract
This work describes the implementation and results related to the numerical error and the temporal analysis of the calibration methodology of an embedded magnetometer. The system tested is part of an attitude determination system for aerospace applications, which uses low-cost tri-axial MEMS (micro electromechanical systems) magnetometer and an ARM Cortex-M4 as the processing unit. In order to validate and develop an efficient methodology for in run calibration, two different calibration algorithms were evaluated: a Two-Step process and a Minimum Variance Estimator (MVE), both described in the literature. The sensor's model consists of nine parameters to be estimated: bias, scale factor and misalignment angles between axes. A Test Unit was developed to evaluate both the numerical and temporal results of calibration. It uses IEEE 754 floating point standard variables with 32-bit and 64-bit wide in the ARM Cortex-M4. The former can be executed within the processor's Float Point Unit and the later cannot. The numerical data are compared to the one obtained using a Matlab code simulation. Result shows that the Two-Step present numerical instability for 32-bit wide variables, while the Minimum Variance Estimator has an inferior temporal execution but did not present such numerical issues.

## Getting Started: Prerequisites

A list containing the softwares and hardware used in this project is available below.

#### Software
* [GNU Arm Embedded Toolchain Version 8-2018-q4-major](https://developer.arm.com/tools-and-software/open-source-software/developer-tools/gnu-toolchain/gnu-rm/downloads)
* [Microsoft Visual Studio Enterprise 2015 Version 14.0.25431.01 Update 3](https://docs.microsoft.com/en-us/visualstudio/releasenotes/vs2015-version-history)
* [Microsoft .NET Framework Version 4.7.03190](https://docs.microsoft.com/en-us/visualstudio/releasenotes/vs2015-version-history)
* [VisualGDB v5.2](https://visualgdb.com/history/)

##### Used Visual GDB Packages
* STM32 Devices v2018.12R2
* OpenOCD v20170609-0.10.0
* Embedded Profiler and Fast Semihosting v3.1
* CppUTest v3.8r1

#### Hardware
[STM32F4DISCOVERY](https://www.st.com/en/evaluation-tools/stm32f4discovery.html)

## Content

This section describes the files shared here.

### Matlab

Inside the "matlab" folder, you will find a script that simulates Figure 4 and the data displayed in Table 5 of the article. To construct the figure, a function that performs a polynomial fit of order n is also available.

### Visual Studio

The most important constributions of this package are the libraries implemented in the header, labeled as "hdr", and source, "src", folders available inside "visual-studio/magnetic-sensor-calib/". Following the alphabetical order, we have:

* magdata.h: the main purpouse behind this library is to load, to perform the tests made, a part of an extesive volume of raw magnetic data from a magnetometer that has gone through a long period inside a Helmholtz coil. For further information about the magnetometer and  Helmholtz coil, see: D. S. Batista, F. Granziera, M. C. Tosin and L. F. de Melo, “[Three-Axial Helmholtz Coil Design and Validation for Aerospace Applications](https://ieeexplore.ieee.org/document/8062793),” in IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 1, pp. 392-403, Feb. 2018.
* matrix.h: implements a bunch of linear algebra functions, e.g., matrix addition, matrix multiplication, and inverse (using the LU decomposition method) of a matrix. 
* mve.h: performs the auto-calibration routine based on the minimum variance estimator (MVE) approach.
* test.h: functions used to run the tests that make use of the CppUTest package.
* twostep.h: performs the auto-calibration routine based on the two-step process approach.


## Running the Tests


## Citing this Repository and License
