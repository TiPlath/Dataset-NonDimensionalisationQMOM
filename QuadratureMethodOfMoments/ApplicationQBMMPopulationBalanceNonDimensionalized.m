%% Quadrature Method of Moments non-dimensionalized Application
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script utilizes the non-dimensionalized quadrature method of moments
% for combined aggregation, breakage and growth mechanisms. The kernels can
% be adjusted if needed It reads in a PSD and evolves its set of moments in
% time. One can compute physical values from these moments to get
% information about the evolved particle size distribution. This script is
% designed to apply the non-dimensionalised model to experimental data by
% fitting the Sauter mean diameter.
% 
% 
% INPUT: A .csv file for the function ReadData, containing information of 
%        a particle size distribution, i.e. diameters (column 1) and
%        probabilities/densities (column 2).
% 
% OUTPUT: iterM     A vector containing all 2N_delta moments for all 
%                   timesteps computed.
%         d32       A vector containing the Sauter mean diameter for all
%                   timesteps

%% MATLAB options
clear all
close all
clc

% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables
% timestep derived from pi_t
dt=0.01;
% max time derived from pi_t
tmax=11;
% Number of time steps
nt=round(tmax/dt) + 1;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,25)
N_delta=3;
% Maximum number of moments
mMax = 2*N_delta;
% store t for every iteration
itert = zeros(nt,1);
% call function to read distribution with units. In our case we read in a
% volume-based number distribution (n_V) [m^-3].
n_V = ReadData('Data/InitialLactoseMCCPVP-n_V.csv',1,0,1);
n_V_N13 = ReadData('Data/N13LactoseMCCPVP-n_V.csv',1,0,1);

% Conversion to get n_L to evolve moments for L to get the right Sauter
% mean diameter in the end
% CDFn_L
CDFn_V = convertPDFtoCDF(n_V);
CDFn_V_N13 = convertPDFtoCDF(n_V_N13);

% Conversion:
k_v = pi/6;
CDFn_L = CDFn_V;
CDFn_L(:,1) = (CDFn_V(:,1)/k_v).^(1/3);
CDFn_L_N13 = CDFn_V_N13;
CDFn_L_N13(:,1) = (CDFn_V_N13(:,1)/k_v).^(1/3);

% Convert CDF to PDF
n_L = convertCDFtoPDF(CDFn_L);
n_L_N13 = convertCDFtoPDF(CDFn_L_N13);
%% Typical Kernel values for non-dimensionalization
% Aggregation rate [1/p*s*m^3]
a_t = zeros(mMax,1);
a_t(:) = 5e4;
% dimensionless kernels based on aggregation timescale (assume that we
% change G_t and beta_t to fit the dimensionless kernel values)
pi_G = zeros(N_delta,1);
pi_G(:) = 0.0152;
pi_beta_t = zeros(N_delta,1);
pi_beta_t(:) = 0.1123;
pi_G_L = zeros(N_delta,1);
pi_G_L(:) = 0.000152;
pi_beta_Lt = zeros(N_delta,1);
pi_beta_Lt(:) = 0.1679;
% Number of fragments from a breakage event
N_f = 8;
%% Define Scale variables (pre-processing) 
% symmetric fragmentation for n_V
b_alpha = @(V_alpha,k,N_f) N_f^(1-k)*V_alpha.^k;
% symmetric fragmentation for n_L
b_alpha_L = @(L_alpha,k,N_f) N_f^((3-k)/3)*L_alpha.^k;
% Compute total number of particles to define particle unit scale P
M_0 =  sum(n_V(1:end-1,2).*diff(n_V(:,1)));
% Compute total volume to define length scale L
M_1 = sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1)));
% Particle scale [p]
U_P = M_0;
% Volume scale [m^3]
U_V = (M_1 * U_P^(-1));
% Time scale [s]
U_a = a_t(1)*U_V;
U_T = U_P^(-1)*U_a^(-1);
pi_a_t = zeros(N_delta,1);
pi_a_t(:) = 1;
%% Define scaled variables (pre-processing)
% define dimensionless particle size V [m^3] and distribution function n_L [p/m^3]
% to ensure we get dimensionless moments
n_V(:,1) = n_V(:,1)./(U_V);
n_V(:,2) = n_V(:,2).*U_V.*U_P^(-1);
n_V_N13(:,1) = n_V_N13(:,1)./(U_V);
n_V_N13(:,2) = n_V_N13(:,2).*U_V.*U_P^(-1);
n_L_N13(:,1) = n_L_N13(:,1)./(U_V^(1/3));
n_L_N13(:,2) = n_L_N13(:,2).*U_V^(1/3).*U_P^(-1);
n_L(:,1) = n_L(:,1)./(U_V^(1/3));
n_L(:,2) = n_L(:,2).*U_V^(1/3).*U_P^(-1);
% dimensionless time
pi_dt=dt;
pi_tmax=tmax;
pi_nt = round(pi_tmax/pi_dt)+1;
pi_t=(0:pi_nt+1);
%Initialize initial values for iteration vectors (useful for debugging and plotting)
% A matrix storing the momenta calculated from n_V for every time vstep
iterM=zeros(nt,2*N_delta);
iterM_L=zeros(nt,2*N_delta);
% store weights and nodes of every iteration
iterV_alpha = zeros(nt,N_delta);
iterW_alpha = zeros(nt,N_delta);
%% get initial weights and nodes
% compute the moments of the initial length based number distribution n_L
M = ComputeMoments(n_V(:,1), n_V(:,2), mMax);
M_N13 = ComputeMoments(n_V_N13(:,1), n_V_N13(:,2), mMax);
M_L_N13 = ComputeMoments(n_L_N13(:,1), n_L_N13(:,2), mMax);
M_L = ComputeMoments(n_L(:,1), n_L(:,2), mMax);
% compute weights and nodes using the Wheeler algorithm
[V_alpha,w_alpha] = Wheeler(M(1:2*N_delta,2*N_delta),N_delta);
[V_alpha_N13,w_alpha_N13] = Wheeler(M_N13(1:2*N_delta,2*N_delta),N_delta);
[L_alpha_N13,wL_alpha_N13] = Wheeler(M_L_N13(1:2*N_delta,2*N_delta),N_delta);
[L_alpha,wL_alpha] = Wheeler(M_L(1:2*N_delta,2*N_delta),N_delta);
%% evolve weights and nodes in time according to QMOM
% Store inital moments, nodes and weights of the volume-based distribution for
% further processing (N_V)
iterM(1,:) = getMomenta(V_alpha,w_alpha);
iterM_N13(1,:) = getMomenta(V_alpha_N13,w_alpha_N13);
iterM_L_N13(1,:) = getMomenta(L_alpha_N13,wL_alpha_N13);
iterM_L(1,:) = getMomenta(L_alpha,wL_alpha);
iterV_alpha(1,1:N_delta) = V_alpha;
iterW_alpha(1,1:N_delta) = w_alpha;
% Time loop, forward Euler
for i=2:nt
    % Compute one time step of QMOM (with N_V)
    [S_V] = GrowthHydrodynamicAggregationPowerLawBreakage(V_alpha,w_alpha,pi_G,pi_a_t,pi_beta_t,b_alpha,N_f);
    % Evolve Moments in time by forward Euler
    iterM(i,:) = iterM(i-1,:)' + S_V.*pi_dt;
    % we compute weights and nodes using the Wheeler algorithm
    [V_alpha,w_alpha] = Wheeler(iterM(i,:)',N_delta);
    % Compute one time step of QMOM (with N_L)
    [S_VL] = GrowthHydrodynamicAggregationPowerLawBreakageLengthBased(L_alpha,wL_alpha,pi_G_L,pi_a_t,pi_beta_Lt,b_alpha_L,N_f);
    % Evolve Moments in time by forward Euler
    iterM_L(i,:) = iterM_L(i-1,:)' + S_VL.*pi_dt;
    % we compute weights and nodes using the Wheeler algorithm
    [L_alpha,wL_alpha] = Wheeler(iterM_L(i,:)',N_delta);
    % Assign weights, nodes and momenta  to iteration vectors
    iterV_alpha(i,1:N_delta) = V_alpha;
    iterW_alpha(i,1:N_delta) = w_alpha;
end
%% Post-Processing
% define stepsize for plotting to have less markers
stepsize = (pi_nt-1)*pi_dt;

%% plot time evolution of all moments
figure(10)
axes10 = axes('Parent', figure(10));
legendLabels = cell(1,size(iterM,2));
for i = 1:mMax
    plot(pi_t(1:stepsize:pi_nt)*dt,iterM(1:stepsize:end,i)/max(iterM(:,i)),'-x', 'MarkerSize', 14,  'LineWidth', 2)
    hold on
    legendLabel{i} = ['$\pi_{\widetilde{M}_{', num2str(i-1, '%i'), '}}$'];
end
xlabel('$\pi_t$','Interpreter','Latex')
ylabel('$\pi_{\widetilde{M}_{k}}$','Interpreter','Latex')
legend(legendLabel{:}, 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes10,'FontSize', 16)

% Reset the color order to repeat the defined colors
set(axes10, 'ColorOrderIndex', 1);
for i = 1:mMax
    plot(pi_t(end-2)*dt,iterM_N13(1,i)/max(iterM_N13(1,:)),'o', 'MarkerSize', 14,  'LineWidth', 2)
end

%% scaling the moments to the same units
for k = 1:mMax-1
    iterM_L(:,k+1) = iterM_L(:,k+1).^(1/(k));
    iterM_L_N13(k+1) = iterM_L_N13(k+1).^(1/(k));
end

for k = 1:mMax-1
    iterM(:,k+1) = iterM(:,k+1).^(1/(3*k));
    iterM_N13(k+1) = iterM_N13(k+1).^(1/(3*k));
end
%% Plot error V
for i = 0:mMax-1
    E_fit(i+1) = abs((iterM((pi_nt),i+1)/iterM_N13(i+1)) - 1);
end

figure(11)
axes11 = axes('Parent', figure(11));
legendLabels = cell(1,size(iterM,2));
stem(E_fit(1,1), "-",'MarkerSize', 14, 'LineWidth',2)
hold on
for i = 3:size(E_fit,2)
    stem(E_fit(1,i), "-",'MarkerSize', 14, 'LineWidth',2)
    hold on
end
xlabel('$\pi_t$','Interpreter','Latex')
ylabel('$\pi_{\widetilde{M}_{k}}$','Interpreter','Latex')
legend('$E_{\widetilde{M}_0}$', '$E_{\widetilde{M}_2}$', '$E_{\widetilde{M}_3}$', '$E_{\widetilde{M}_4}$', '$E_{\widetilde{M}_5}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes11,'FontSize', 16)
fprintf('relative mean moment error for the V10 fit is: E = %g\n', sum(E_fit)/length(E_fit))
%% Plot error L
for i = 0:mMax-1
    E_L_fit(i+1) = abs((iterM_L((pi_nt),i+1)/iterM_L_N13(i+1)) - 1);
end

figure(12)
axes12 = axes('Parent', figure(12));
legendLabels = cell(1,size(iterM,2));
stem(E_L_fit(1,1), "-",'MarkerSize', 14, 'LineWidth',2)
hold on
for i = 3:size(E_L_fit,2)
    stem(E_L_fit(1,i), "-",'MarkerSize', 14, 'LineWidth',2)
    hold on
end
xlabel('$\pi_t$','Interpreter','Latex')
ylabel('$\pi_{\widetilde{M}_{k,L}}$','Interpreter','Latex')
legend('$E_{\widetilde{M}_{0,L}}$', '$E_{\widetilde{M}_{2,L}}$', '$E_{\widetilde{M}_{3,L}}$', '$E_{\widetilde{M}_{4,L}}$', '$E_{\widetilde{M}_{5,L}}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes12,'FontSize', 16)
fprintf('relative mean moment error for the d32 fit is: E = %g\n', sum(E_L_fit)/length(E_L_fit))
%% plot d32 fit

d32_N13 = (iterM_L_N13(4))/(iterM_L_N13(3)) * U_V^(1/3);
d32 = zeros(1,pi_nt);
d32 = (iterM_L(:,4))./(iterM_L(:,3)) * U_V^(1/3);

figure(13)
axes13 = axes('Parent', figure(13));
plot(pi_t(1:stepsize:pi_nt)*dt,d32(1:stepsize:pi_nt), "x", 'MarkerSize', 14, 'LineWidth',2)
hold on
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$d_{32}$ [m]','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes13,'FontSize', 16)
line([axes13.XLim(1), axes13.XLim(2)], [d32_N13, d32_N13], 'Color', 'red', 'LineStyle', '-', 'LineWidth', 2);
legend('$d_{32}$', '$d_{32}^{N13}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
%% plot V10

V10_N13 = (iterM_N13(2))/(iterM_N13(1)) * U_V^(3/3);
V10 = zeros(1,pi_nt);
V10 = (iterM(:,2))./(iterM(:,1)) * U_V^(3/3);

figure(14)
axes14 = axes('Parent', figure(14));
plot(pi_t(1:stepsize:pi_nt)*dt*U_T,V10(1:stepsize:pi_nt), "x", 'MarkerSize', 14, 'LineWidth',2)
hold on
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$V_{10}$ [m\textsuperscript{3}]','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes14,'FontSize', 16)
line([axes14.XLim(1), axes14.XLim(2)], [V10_N13, V10_N13], 'Color', 'red', 'LineStyle', '-', 'LineWidth', 2);
legend('$V_{10}$', '$V_{10}^{N13}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')

%% plot moments against moment order semilogarithmic V
momentOrder = 0:1:(size(iterM,2) - 1);
figure(15)
axes15 = axes('Parent', figure(15));
plot(momentOrder,iterM(end,:), '--x','MarkerSize', 14, 'LineWidth',2)
hold on
plot(momentOrder,iterM_N13(:),'--x','MarkerSize', 14, 'LineWidth',2)
xlabel('Moment order','Interpreter','Latex')
ylabel('$\pi_{\bar{M}_{k}}$','Interpreter','Latex')
legend('$\pi_{\bar{M}_k}(\pi_t=11)$', '$\pi_{\bar{M}_k} (N13)$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes15,'FontSize', 16)
%% plot moments against moment order semilogarithmic L

momentOrder = 0:1:(size(iterM,2) - 1);
figure(16)
axes16 = axes('Parent', figure(16));
plot(momentOrder,iterM_L(end,:), '--x','MarkerSize', 14, 'LineWidth',2)
hold on
semilogy(momentOrder,iterM_L_N13(:),'--x','MarkerSize', 14, 'LineWidth',2)
xlabel('Moment order','Interpreter','Latex')
ylabel('$\pi_{\bar{M}_{k,L}}$','Interpreter','Latex')
legend('$\pi_{\bar{M}_{k,L}}(\pi_t=11)$', '$\pi_{\bar{M}_{k,L}} (N13)$', 'location','northwest', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes16,'FontSize', 16)