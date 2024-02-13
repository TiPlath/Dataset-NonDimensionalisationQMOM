%% Converts length based cumulative volume distribution to a volume-based
%% number density distribution (n_V)
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script is used to cut particle size distributions with high size
% ratios to prevent inaccuracies. It cuts a certain quantile at the minimum
% and maximum particle size.
% 
% INPUT:    CDF     A cumulative probability density distribution function, containing 
%                   i.e. diameters (column 1) and probabilities/densities (column 2).
% 
% OUTPUT:   cutCDF  The cut particle size distribution


function cutCDF = cutHighSizeRatio(CDF,quantileMin,quantileMax,minPolydispersity)

if quantileMin == 0 && quantileMax == 1
    cutCDF = CDF;
    return;
end
% maximum probability
maxProb = CDF(end,2);
% normalize
CDF(:,2) = CDF(:,2)/max(CDF(:,2));
% Size ratio = rMax/rMin
SR = CDF(end,1)/CDF(1,1);
sprintf('Initial size ratio is %f', SR)
radiusMin = getRadiusByQuantile(CDF,quantileMin);
radiusMax = getRadiusByQuantile(CDF,quantileMax);
% minimum polyidispersity at the base
radiusMinCut = min(radiusMin*(1+minPolydispersity),radiusMax);
% cut off min
while CDF(1,1) <= radiusMinCut
    CDF(1,:) = [];
end
CDF = [radiusMinCut,quantileMin;CDF];
if quantileMin ~= 0
    CDF = [radiusMin,0;CDF];
end
% cut off max
while CDF(end,1) >= radiusMax
    CDF(end,:) = [];
end
radiusMaxCut = max(radiusMax-(1-minPolydispersity)*(radiusMax-CDF(end,1)),radiusMin);

% delete columns above probability 1.0
while CDF(end,2) >= 1.0
    CDF(end,:) = [];
end
% delete end if radius above radiusMax
if CDF(end,1) >= radiusMax
    CDF(end,:) = [];
end

SR = CDF(end,1)/CDF(1,1);
sprintf('Cut size ratio is %f', SR)
%de-normalize to not loose information
CDF(:,2) = CDF(:,2)*maxProb;
cutCDF = CDF;
end

function radius = getRadiusByQuantile(PDF,quantile)
if quantile > 1 || quantile < 0
    error('quantile has to be between 0 and 1')
end
%Get readius index nearest to quantile
[~,nearestIndex] = min(abs(PDF(:,2)-quantile));
radius = PDF(nearestIndex,1);
end