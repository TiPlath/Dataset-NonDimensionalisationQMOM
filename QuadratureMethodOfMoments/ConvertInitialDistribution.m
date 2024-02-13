%% Converts length based cumulative volume distribution to a volume-based
%% number density distribution (n_V)
% Author = Plath, Timo & Weinhart, Thomas
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script is used to convert the length-based particle size
% distribution. It reads in the lenght-based particle size distribution and
% converts it to a volume-based particle size distribution. At the end of
% every conversion the distribution is plotted and the volume is checked
% for conservation.
% 
% INPUT:    A .csv file for the function ReadData, containing information 
%           of a particle size distribution, i.e. diameters (column 1) and
%           probabilities/densities (column 2).
% 
% OUTPUT:   A plot of the inital volume-based particle size distribution.
%           And a matrix which one can extract into a .csv file to get the
%           volume-based number particle size distribution for the
%           simulations.
%% MATLAB options
clear all
close all
clc

% The volume extracted from a QICPIC q3[mm^-1] measurement for the initial
% Lactose, MCC, PVP mixture of the dataset from Plath et al.
V = 2.05e-6;
fprintf('bulk volume based on q3 distribution from Plath et al.: %g m^-3\n', V)
% CVDF_L
% Interpretation: CVDF_L(i,2) is the bulk volume of particles of
% diameter less than CVDF_L(i,1); normalised to 1m^3

% Read Q3:
CDFv_L = ReadData('Data/InitialLactoseMCCPVP-CDF-v_L.csv',1e6,0,0.995);
% normalize data to 1. The initial data is already normalized by volume to 100%
CDFv_L(:,2) = CDFv_L(:,2)/max(CDFv_L(:,2));
% Read Q3 of N13:
CDFv_L_N13 = ReadData('Data/N13LactoseMCCPVP-CDF-v_L.csv',1e6,0,1);
% normalize data to 1. The initial data is already normalized by volume to 100%
CDFv_L_N13(:,2) = CDFv_L_N13(:,2)/max(CDFv_L_N13(:,2));

% Plot Q3 before granulation with Q3 of N13 after granulation
figure(1)
axes1 = axes('Parent', figure(1));
plot(CDFv_L(:,1),CDFv_L(:,2),'-x', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot(CDFv_L_N13(:,1),CDFv_L_N13(:,2), '-x', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('$L$ [m]','Interpreter','Latex')
ylabel('$Q_3$ [m\textsuperscript{2}]','Interpreter','Latex')
legend( '$Q_3^{initial}$','$Q_3^{N13}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes1,'FontSize', 16)

% Plot: 
figure(2)
subplot(2,3,1)
semilogx(CDFv_L(:,1),CDFv_L(:,2))
xlabel('L [m]'); ylabel('CDFv_L [m^2]'); 

% Check:
fprintf('Normalised bulk volume based on CDFv_L: %g m^3\n', CDFv_L(end,2))

% CDFv_V
% Interpretation: CDFv_V(i,2) is the bulk volume of particles of
% volume less than CDFv_V(i,1) 

% Conversion:
k_v = pi/6;
CDFv_V = CDFv_L;
CDFv_V(:,1) = k_v .* CDFv_L(:,1).^3;

% Plot: 
subplot(2,3,2)
semilogx(CDFv_V(:,1),CDFv_V(:,2))
xlabel('V [m^3]'); ylabel('CDFv_V [m^3]'); 

% Check:
fprintf('Normalised bulk volume based on CDFv_V: %g m^3\n', CDFv_V(end,2))

% v_V
% v_V(i,2)*(v_V(i+1,1)-v_V(i,1)) is the bulk volume fraction of
% particles of volume between v_V(i,1) and v_V(i+1,1)

% Conversion:
v_V = CDFv_V;
v_V(:,2) = [diff(CDFv_V(:,2))./diff(v_V(:,1)); 0];

% Plot: 
subplot(2,3,4)
stairs(v_V(:,1),v_V(:,2))
set(gca,'XScale','log'); xlabel('V [m^3]'); ylabel('v_V [m^3/m^3]'); 

% Check: 0th moment == volume
fprintf('Normalised bulk volume based on v_V: %g m^3\n', sum(v_V(1:end-1,2).*diff(v_V(:,1))))

% v_L
% Conversion:
v_L = CDFv_L;
v_L(:,2) = [diff(CDFv_L(:,2))./diff(v_L(:,1)); 0];

% Check: 0th moment == volume
fprintf('Normalised bulk volume based on v_L: %g m^3\n', sum(v_V(1:end-1,2).*diff(v_V(:,1))))

% n_L
% Conversion:
n_L = v_L;
n_L(:,2) = v_L(:,2)./v_L(:,1).^3;

% Plot:
subplot(2,3,6)
stairs(n_L(1:end-1,1),n_L(1:end-1,2))
set(gca,'XScale','log'); set(gca,'YScale','log'); xlabel('L [m]'); ylabel('n_L [P/m]'); 


% Check: 3rd moment == volume
fprintf('Normalised bulk volume based on n_L: %g m^3\n', sum(n_L(1:end-1,1).^3.*n_L(1:end-1,2).*diff(n_L(:,1))))

% n_V
% n_V(i,2)*(n_V(i+1,1)-n_V(i,1)) is the number fraction of
% particles of volume between PNDF_V(i,1) and PNDF_V(i+1,1)
% Conversion: Here we have to do assumptions on how the continuous PDF's
% look; if we assume all particles have the smallest possible radius, then 
% int(v*n_v_i,v_i,v_i+1) = int(v_v,v_i,v_i+1) 
% int(v_i*n_v_i,v_i,v_i+1) = int(v_v,v_i,v_i+1) 
% n_v_i*v_i*(v_i+1-v_i) = v_v_i*(v_i+1-v_i) 
% n_v_i = 2*v_v_i/v_i
% note: a better assumption would be that the n_V is piecewise linear,
% but then the conversion formula is more complex
n_V = v_V;
n_V(:,2) = v_V(:,2)./v_V(:,1);

% Plot: 
subplot(2,3,5)
stairs(n_V(1:end-1,1),n_V(1:end-1,2))
set(gca,'XScale','log'); set(gca,'YScale','log'); xlabel('V [m^3]'); ylabel('n_V [P/m^3]'); 

% Check: 1st moment == volume
fprintf('Normalised bulk volume based on n_V: %g m^3\n', sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1))))


% re-normalize
n_V(:,2) = n_V(:,2)*V;
% 1st moment of re-normalized PSD == Volume
fprintf('Total particle volume based on re-normalized n_V: %g m^3\n', sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1))))
% 0th moment of re-normalized PSD == number
fprintf('Total particle number based on re-normalized n_V: %g P\n', sum(n_V(1:end-1,2).*diff(n_V(:,1))))

fprintf('Conversion check to get the right Sauter mean diameter (d_32):\n')
% CDFn_L
CDFn_V = convertPDFtoCDF(n_V);

% Check:
fprintf('Total particle number based on CDFn_V: %g P\n', CDFn_V(end,2))

% Conversion:
CDFn_L = CDFn_V;
CDFn_L(:,1) = CDFn_V(:,1).^(1/3);

% Plot: 
subplot(2,3,3)
semilogx(CDFn_L(:,1),CDFn_L(:,2)/CDFn_L(end,2))
xlabel('L [m]'); ylabel('CDFn_L [m^3]');

% Check:
fprintf('Total particle number based on CDFn_L: %g P\n', CDFn_L(end,2))


% n_L 
% Conversion
n_L = convertCDFtoPDF(CDFn_L);

% Check: 1st moment == particle number
fprintf('Total particle number based on n_L: %g P\n', sum(n_L(1:end-1,1).^0.*n_L(1:end-1,2).*diff(n_L(:,1))))
% Check: 3rd moment == volume
fprintf('Total bulk volume based on n_L: %g m^3\n', sum(n_L(1:end-1,1).^3.*n_L(1:end-1,2).*diff(n_L(:,1))))


