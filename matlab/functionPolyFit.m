function [fittpar,Rsq,Rsq_adj] = functionPolyFit(x,y,n) 
%Perform the n-order fit of y and x, returning the fitting parameters and 
%coefficients of determination (normal and ajusted ones).
%
%This Matlab function is used in the article:
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
%@Inputs:
%   x: Independent variable lying on the abscissa axis.
%   y: Dependent variable lying on the ordinate axis.
%	n: Order of the polynomial fit.
%
%@Outputs:
%	fittpar: Fitting parameters.
%	Rsq: Coefficient of determination.
%	Rsq_adj: Adjusted coefficient of determination.
%

    %Compute the paramenters from a fit of order n 
    fittpar = polyfit(x,y,n);
    
    %Compute the dependent variable obtained from the fitted equation
    yfit = polyval(fittpar,x);
    
    %Compute the residual
    yresid = y-yfit;
    
    %Compute the sum of squared (SS) errors and also its total value based
    %on y true 
    SSresid = sum(yresid.^2);
    SStotal = var(y)*(length(y)-1);
    
    %Compute the coefficients of determination
    Rsq = 1-(SSresid/SStotal);
    Rsq_adj = 1-((SSresid/SStotal)*(length(y)-1)/(length(y)-length(fittpar))); % correct the coefficient based on the true value
    
end
