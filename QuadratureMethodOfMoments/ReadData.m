%% Pre-processing of data
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Reads in particle size distributions from data files and scales them to
% their respective units by a sizeScaling factor and cuts them if necessary
% 
% INPUT:  fileName      strings of the .csv file
%         sizeScaling   double for conversion into SI Units
%         QuantileMin   double determining the quantile to cut the left 
%                       hand side of a distribution
%         QuantileMax   double determining the quantile to cut the right
%                       hand side of a distribution
% 
% OUTPUT: T             A vector containing the diameters and 
%                       probabilities/densities of the .csv file

function T = ReadData(fileName,sizeScaling,QuantileMin,QuantileMax)
% Define scaling for conversion to SI-Units (size in [m] and probabilities
% in [0,1] [%] where the type of distribution determines the weight of
% percentages
% read initial Lactose distribution
T = readtable(fileName,'PreserveVariableNames',true, 'Delimiter', ',');
% convert table into an array to make it more readable
T = table2array(T);
% convert into SI Units
T(:,1) = T(:,1)./sizeScaling;
% Cut the PSD at a certain quantile to avoid inaccuracies due to high size
% ratios
T = cutHighSizeRatio(T,QuantileMin,QuantileMax,0.1);
end