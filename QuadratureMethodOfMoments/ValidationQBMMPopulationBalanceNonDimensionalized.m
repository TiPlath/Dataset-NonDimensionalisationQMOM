%% Quadrature Method of Moments non-dimensionalized Validation
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script utilizes the non-dimensionalized quadrature method of moments
% for combined aggregation, breakage and growth mechanisms. The kernels can
% be adjusted if needed It reads in a PSD and evolves its set of moments in
% time. One can compute physical values from these moments to get
% information about the evolved particle size distribution. This script is
% specifically designed to validate the non-dimensionalised model for pure
% growth, pure aggregation as well as aggregation and breakage. Therefore
% the user has to set pi_beta and pi_G to specific values.
% 
% 
% INPUT: A .csv file for the function ReadData, containing information of 
%        a particle size distribution, i.e. diameters (column 1) and
%        probabilities/densities (column 2).
% 
% OUTPUT: iterM     A vector containing all 2N_delta moments for all 
%                   timesteps computed.

%% MATLAB options
clear all
close all
clc

% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables
combinations = 5;
% timestep derived from pi_t
dt=0.1;
% max time derived from pi_t
tmax=10;
% Number of time steps
nt=round(tmax/dt) + 1;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,25)
N_delta=1;
% Maximum number of moments
mMax = 2*N_delta;
% store M0,M1,t for every iteration with different aggregation rates
iterM0 = zeros(nt,1);
iterM1 = zeros(nt,1);
itert = zeros(nt,1);
% call function to read distribution with units. In our case we read in a
% volume-based number distribution (n_V) [m^-3].
n_V = ReadData('Data/InitialLactoseMCCPVP-n_V.csv',1,0,1);
%% Typical Kernel values for non-dimensionalization
% Aggregation rate [p/s]
a = zeros(mMax,1);
a(:) = 1e-8;
% dimensionless kernels based on aggregation timescale (assume that we
% change g_0 and beta_0 to fit the dimensionless kernel values)
pi_G = zeros(N_delta,1);
pi_G(:) = 1;
pi_beta = zeros(N_delta,1);
pi_beta(:) = 0.5;
%% Define Scale variables (pre-processing)
% symmetric fragmentation for n_V
b_alpha = @(V_alpha,k) 2^(1-k)*V_alpha.^k;
% Compute total number of particles to define particle unit scale P
M_0 =  sum(n_V(1:end-1,2).*diff(n_V(:,1)));
% Compute total volume to define length scale L
M_1 = sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1)));
% Particle scale [p]
U_P = M_0;
% Time scale [s]
U_T = a(1)^(-1) * U_P^(-1);
pi_a = zeros(N_delta,1);
pi_a(:) = 1;
% Volume scale [m^3]
U_V = (M_1 * U_P^(-1));
%% Define scaled variables (pre-processing)
% define dimensionless particle size l [m] and distribution function n_L [p/m]
% to ensure we get dimensionless moments
n_V(:,1) = n_V(:,1)./(U_V);
n_V(:,2) = n_V(:,2).*U_V.*U_P^(-1);
% dimensionless time
pi_dt=dt;
pi_tmax=tmax;
pi_nt = round(pi_tmax/pi_dt)+1;
pi_t=(0:pi_nt+1);
%Initialize initial values for iteration vectors (useful for debugging and plotting)
% A matrix storing the momenta calculated from n_V for every time vstep
iterM=zeros(nt,2*N_delta);
% store weights and nodes of every iteration
iterV_alpha = zeros(nt,N_delta);
iterW_alpha = zeros(nt,N_delta);
%% get initial weights and nodes
% compute the moments of the initial length based number distribution n_L
M = ComputeMoments(n_V(:,1), n_V(:,2), mMax);
% compute weights and nodes using the Wheeler algorithm
[V_alpha,w_alpha] = Wheeler(M(1:2*N_delta,2*N_delta),N_delta);
%% evolve weights and nodes in time according to QMOM
% Store inital moments, nodes and weights of the volume-based distribution for
% further processing (N_V)
iterM(1,:) = getMomenta(V_alpha,w_alpha);
iterV_alpha(1,1:N_delta) = V_alpha;
iterW_alpha(1,1:N_delta) = w_alpha;
iterM0(1,1) = iterM(1,1);
iterM1(1,1) = iterM(1,2);
% A matrix to plot the time for different iterations
itert = pi_t(1:pi_nt);
% Time loop, forward Euler
for i=2:nt
    % Compute one time step of QMOM (with N_V)
    [S_V] = GrowthAggregationBreakage(V_alpha,w_alpha,pi_G,pi_a,pi_beta,b_alpha);
    % Evolve Moments in time by forward Euler
    iterM(i,:) = iterM(i-1,:)' + S_V.*pi_dt;
    % we compute weights and nodes using the Wheeler algorithm
    [V_alpha,w_alpha] = Wheeler(iterM(i,:)',N_delta);
    % Assign weights, nodes and momenta  to iteration vectors
    iterV_alpha(i,1:N_delta) = V_alpha;
    iterW_alpha(i,1:N_delta) = w_alpha;
    iterM0(i) = iterM(i,1);
    iterM1(i) = iterM(i,2);
end
%% Post-Processing
% define stepsize for plotting to have less markers
stepsize = (pi_nt-1)*pi_dt;
%% Analytical solutions for Aggregation and breakage pi_M_0, pi_M_1
syms M0(t_) M1(t_)
M0(t_) = (2*pi_beta(1)*iterM(1,1)*exp(pi_beta(1)*t_))/(2*pi_beta(1)-pi_a*iterM(1,1)+iterM(1,1)*pi_a*exp(pi_beta(1)*t_));
M1(t_) = iterM(1,2) - 2*pi_G(1)*log(2*pi_beta(1))/pi_a + 2*pi_G(1)*log(2*pi_beta(1)-pi_a(1)*iterM(1,1)+iterM(1,1)*pi_a*exp(pi_beta(1)*t_))/pi_a;

if pi_beta(1) == 0
    M0(t_) = 2*iterM(1,1)/(2+pi_a(1)*iterM(1,1)*t_);
end
if pi_G == 0
    M1(t_) = iterM(1,2);
end
%% plot time evolution of all moments
figure(10)
axes10 = axes('Parent', figure(10));
legendLabels = cell(1,size(iterM,2));
for i = 1:2
    plot(pi_t(1:stepsize:pi_nt)*dt,iterM(1:stepsize:end,i),'x', 'MarkerSize', 14,  'LineWidth', 2)
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

fplot(subs(M0),[0 pi_tmax], 'LineWidth', 2, 'DisplayName', '$\pi_{\widetilde{M}_0}^{analytic}$');
hold on
fplot(subs(M1),[0 pi_tmax], 'LineWidth', 2, 'DisplayName', '$\pi_{\widetilde{M}_1}^{analytic}$')

if pi_G(1) == 0 && pi_beta(1) == 0.25
    ylim([0 1])
end
