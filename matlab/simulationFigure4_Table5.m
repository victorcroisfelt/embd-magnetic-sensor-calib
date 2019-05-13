%This Matlab script can be used to generate all figures in the article:
%
%Victor Croisfelt Rodrigues, Guilherme Brandão da Silva, Daniel Strufaldi
%Batista, Marcelo Carvalho Tosin, and Francisco Granziera Junior, 
%"Numerical Error and Temporal Analysis of Embedded Magnetic Sensor 
%Calibration using ARM Cortex-M4 Processor," 2018 13th IEEE International
%Conference on Industry Applications (INDUSCON), São Paulo, Brazil, 2018,
%pp. 1387-1394.
%
%Download paper: https://doi.org/10.1109/INDUSCON.2018.8627205
%
%This is version 2.0 (Last edited: 04-15-2019)
%
%License: This code is licensed under the GPLv3 license. If you in any way
%use this code for research that results in publications, please reference 
%our original article as shown above.
%

%Initialization
close all;
clearvars;

%% STM32407VG specifications

%Define the clock frequency (Hz) of the board 
boardClkFreq = 168e6;

%Obtain the period (s) 
boardPeriod = 1/boardClkFreq;

%Define the evaluated range of the number of samples used to obtain the 
%calibration parameters
nbrOfEvaluatedSamples = 100:25:325;

%% Two-step process: assigning collected data

%Legend:
%   SP: single-precision number (FPU On)
%   DP: double-precision number
%

%Save the adquired number of clock cycles (cc) from different experiments realized on the board
TS.SP.cc = [564356 682812 811821 938876 1065811 1192703 1320155 1442925 1569526 1696780];
TS.SP.p = 1e3*TS.SP.cc*boardPeriod; % period (p) in ms

TS.DP.cc = [4261460 5243867 6234717 7223712 8208492 9197913 10202235 11157172 12173356 13126484];
TS.DP.p = 1e3*TS.DP.cc*boardPeriod; % period (p) in ms

%% Minimum variance estimator: assigning collected data

%Legend:
%   SP: single-precision number (FPU On)
%   DP: double-precision number
%

%Save the adquired number of clock cycles (cc) from different experiments realized on the board
MVE.SP.cc = [870018 1071007 1270468 1471545 1667509 1868669 2069493 2267170 2465332 2663186];
MVE.SP.p= 1e3*MVE.SP.cc*boardPeriod; % period (p) in ms

MVE.DP.cc = [4795924 5919639 7000392 8106636 9208770 10309309 11426057 12555940 13655069 14743104];
MVE.DP.p = 1e3*MVE.DP.cc*boardPeriod; % period (p) in ms

%% Plotting Figure 4

figure
hold on; box on;

plot(nbrOfEvaluatedSamples,TS.SP.p,'d--');
plot(nbrOfEvaluatedSamples,MVE.SP.p,'^--');
plot(nbrOfEvaluatedSamples,TS.DP.p,'o--');
plot(nbrOfEvaluatedSamples,MVE.DP.p,'*--');
    
ylim([0 100])
xlim([100 325]);
xticks(100:25:325);

ylabel('Execution time (ms)')
xlabel('Samples number')  

legend('IEEE 32-bits: Two-Step (FPU On)','IEEE 32-bits: MVE (FPU On)','IEEE 64-bits: Two-Step','IEEE 64-bits: MVE','location','NorthWest');
    
%% Fitting data: Table V

%Legend:
%   fp: fitting parameters
%   coeffDet: coefficient of determination (R^2)
%   adjCoeffDet: adjusted coefficient of determination (R^2)adj
%

%Two step process
[TS.SP.fp,TS.SP.coeffDet,TS.SP.adjCoeffDet] = functionPolyFit(nbrOfEvaluatedSamples,TS.SP.p,1); % single precision
[TS.DP.fp,TS.DP.coeffDet,TS.DP.adjCoeffDet] = functionPolyFit(nbrOfEvaluatedSamples,TS.DP.p,1); % double precision

%MVE
[MVE.SP.fp,MVE.SP.coeffDet,MVE.SP.adjCoeffDet] = functionPolyFit(nbrOfEvaluatedSamples,MVE.SP.p,1); % single precision
[MVE.DP.fp,MVE.DP.coeffDet,MVE.DP.adjCoeffDet] = functionPolyFit(nbrOfEvaluatedSamples,MVE.DP.p,1); % double precision

%Table V preview
disp('Table V')

disp(['Two-Step ',num2str(TS.SP.fp(1),'%-5.4f'),'xN+',num2str(TS.SP.fp(2),'%-5.4f'),' ',num2str(TS.SP.adjCoeffDet,'%-5.4f'),'  ',num2str(TS.DP.fp(1),'%-5.4f'),'xN+',num2str(TS.DP.fp(2),'%-5.4f'),' ',num2str(TS.DP.adjCoeffDet,'%-5.4f')])
disp(['MVE      ',num2str(MVE.SP.fp(1),'%-5.4f'),'xN+',num2str(MVE.SP.fp(2),'%-5.4f'),' ',num2str(MVE.SP.adjCoeffDet,'%-5.4f'),'  ',num2str(MVE.DP.fp(1),'%-5.4f'),'xN+',num2str(MVE.DP.fp(2),'%-5.4f'),' ',num2str(MVE.DP.adjCoeffDet,'%-5.4f')])
