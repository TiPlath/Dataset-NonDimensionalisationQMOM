%% Quadrature Method of Moments non-dimensionalised
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script utilizes the non-dimensionalized quadrature method of moments
% for combined aggregation, breakage and growth mechanisms. We compute the 
% dimensionalised values on scaled time scales. The kernels can be adjusted 
% if needed It reads in a PSD and evolves its set of moments in time. One
% can compute physical values from these moments to get information about 
% the evolved particle size distribution.
% 
% 
% INPUT: A .csv file for the function ReadData, containing information of 
%        a particle size distribution, i.e. diameters (column 1) and
%        probabilities/densities (column 2).
% 
% OUTPUT: iterM0    A vector containing the zeroth moment for all
%                   combinations and timesteps computed.
%         iterM1    A vector containing the first moment for all
%                   combinations and timesteps computed.
%         iterpi_M  A vector containing all 2*N_delta moment for all
%                   timesteps computed.

%% MATLAB options
clear all
close all
clc

% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables
combinations = 5;
% dimensionless time
pi_dt=0.1;
pi_tmax=10;
pi_nt = round(pi_tmax/pi_dt)+1;
pi_t=(0:pi_nt);
% Number of dirac-delta distributed classes (weights and nodes, 1,...,25)
N_delta=1;
% Maximum number of moments
mMax = 2*N_delta;
% store M0,M1,t for every iteration with different aggregation rates
iterM0 = zeros(pi_nt,combinations);
iterM1 = zeros(pi_nt,combinations);
itert = zeros(pi_nt,combinations);
%% dimensional solution of Eq. 8 for all combinations
% start looping over different kernel value combinations
for j = 1:combinations
    % call function to read distribution with units. In our case we read in a
    % volume-based number distribution (n_V) [m^-3].
    n_V = ReadData('Data/InitialLactoseMCCPVP-n_V.csv',1,0,1);
    %% Typical Kernel values for non-dimensionalization
    % Growth rate [m^3/s]
    G = ones(N_delta,combinations);
    for l = 1:size(G)
        G(l,:) = [1e-14, 1e-15, 1e-16, 1e-17, 1e-18];
    end
    % Aggregation rate [p/s]
    a = zeros(N_delta,combinations);
    for l = 1:size(a)
        a(l,:) = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12];
    end
    % Breakage rate[1/s]
    beta = zeros(N_delta,combinations);
    % CHANGE VALUES??
    for l = 1:size(beta)
        beta(l,:) = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
    end
    % symmetric fragmentation for n_V
    b_alpha = @(V_alpha,k) 2^(1-k)*V_alpha.^k;
    %% define different scales and scale the time
    % Compute total number of particles to define particle unit scale P
    M_0 =  sum(n_V(1:end-1,2).*diff(n_V(:,1)));
    % Compute total volume to define length scale L
    M_1 = sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1)));
    % Particle scale [p]
    U_P = M_0;
    % Time scale [s]
    U_t = a(1,j)^(-1) * U_P^(-1);
    pi_a = zeros(N_delta,1);
    pi_a(:,1) = 1;
    % Volume scale [m^3]
    U_V = (M_1 * U_P^(-1));
    % timestep derived from pi_t
    dt=pi_dt*U_t;
    % max time derived from pi_t
    tmax=pi_tmax*U_t;
    % Number of time steps
    nt=round(tmax/dt)+1;
    % time
    t=(0:nt);
    %% get initial weights and nodes
    % compute the moments of the initial length based number distribution n_L
    M = ComputeMoments(n_V(:,1), n_V(:,2), mMax);
    % compute weights and nodes using the Wheeler algorithm
    [V_alpha,w_alpha] = Wheeler(M(1:2*N_delta,2*N_delta),N_delta);
    %% evolve abscissae and xweights in time according to DQMOM
    % Store inital moments, time, nodes and weights of the Volume based distribution
    % for further processing (N_V)
    iterM(1,:) = getMomenta(V_alpha,w_alpha);
    iterV_alpha(1,1:N_delta) = V_alpha;
    iterW_alpha(1,1:N_delta) = w_alpha;
    iterM0(1,j) = iterM(1,1);
    iterM1(1,j) = iterM(1,2);
    itert(:,j) = t(1:nt);
    % Time loop, forward Euler
    for i=2:nt
        % Compute one time step of DQMOM (with N_V)
        [S_V] = GrowthAggregationBreakage(V_alpha,w_alpha,G(:,j),a(:,j),beta(:,j),b_alpha);
        % Evolve Moments in time by forward Euler
        iterM(i,:) = iterM(i-1,:)' + S_V.*dt;
        % we compute weights and nodes using the Wheeler algorithm
        [V_alpha,w_alpha] = Wheeler(iterM(i,:)',N_delta);
        % Assign weights, nodes and momenta to iteration vectors
        iterV_alpha(i,1:N_delta) = V_alpha;
        iterW_alpha(i,1:N_delta) = w_alpha;
        iterM0(i,j) = iterM(i,1);
        iterM1(i,j) = iterM(i,2);
    end
end
%% non-dimensional solution of Eq. 11
% Define scaled variables (pre-processing)
% define dimensionless particle size l [m] and distribution function n_L [p/m]
% to ensure we get dimensionless moments
n_V(:,1) = n_V(:,1)./(U_V);
n_V(:,2) = n_V(:,2).*U_V.*U_P^(-1);
% dimensionless kernels based on aggregation timescale
pi_G = zeros(N_delta,1);
pi_G(:,1) = G(:,j)*U_V^(-1)*U_t;
pi_beta = zeros(N_delta,1);
pi_beta(:,1) = beta(:,j)*U_t;
%Initialize initial values for iteration vectors (useful for debugging and plotting)
% A matrix storing the momenta calculated from n_V for every time vstep
iterpi_M=zeros(nt,2*N_delta);
% store weights and nodes of every iteration
iterpi_V_alpha = zeros(nt,N_delta);
iterpi_W_alpha = zeros(nt,N_delta);
%% get initial weights and nodes
% compute the moments of the initial length based number distribution n_L
M = ComputeMoments(n_V(:,1), n_V(:,2), mMax);
% compute weights and nodes using the Wheeler algorithm
[V_alpha,w_alpha] = Wheeler(M(1:2*N_delta,2*N_delta),N_delta);
%% evolve weights and nodes in time according to QMOM
% Store inital moments, nodes and weights of the volume-based distribution for
% further processing (n_V)
iterpi_M(1,:) = getMomenta(V_alpha,w_alpha);
iterpi_V_alpha(1,1:N_delta) = V_alpha;
iterpi_W_alpha(1,1:N_delta) = w_alpha;
iterpi_t(:) = pi_t(1:nt);
% Time loop, forward Euler
for i=2:nt
    % Compute one time step of QMOM (with n_V)
    [S_V] = GrowthAggregationBreakage(V_alpha,w_alpha,pi_G,pi_a,pi_beta,b_alpha);
    % Evolve Moments in time by forward Euler
    iterpi_M(i,:) = iterpi_M(i-1,:)' + S_V.*pi_dt;
    % we compute weights and nodes using the Wheeler algorithm
    [V_alpha,w_alpha] = Wheeler(iterpi_M(i,:)',N_delta);
    % Assign weights, nodes and momenta  to iteration vectors
    iterV_alpha(i,1:N_delta) = V_alpha;
    iterW_alpha(i,1:N_delta) = w_alpha;
end
%% Post-Processing
% convert moments to non-dimensional moments
iterpi_M0 = iterM0(:,:).*U_P^(-1);
iterpi_M1 = iterM1(:,:).*U_V^(-1)*U_P^(-1);
% define stepsize for plotting to have less markers
stepsize = pi_nt*pi_dt;
% plot moment pi_M_0 against pi_t for all combinations
marker = ['x','o','+','d','v'];
figure(7)
axes7 = axes('Parent', figure(7));
hold on
for j = 1:5
    tPlot = itert(:,j);    
    iterM0Plot = iterpi_M0(1:stepsize:size(iterpi_M0,1),j);
    plot(tPlot(1:stepsize:size(tPlot,1)) * pi_dt,iterM0Plot, marker(j), 'MarkerSize', 12,  'LineWidth', 2)
end
tPlot = iterpi_t(:);
pi_M1Plot = iterpi_M(1:stepsize:size(iterpi_M,1),1);
plot(tPlot(1:stepsize:size(tPlot,1)) * pi_dt,pi_M1Plot, 's', 'MarkerSize', 12, 'LineWidth', 2)
% reset color order to plot last value
set(axes7, 'ColorOrderIndex', 1);
for j = 1:5
    plot(tPlot(end) * pi_dt, iterpi_M0(end,j), marker(j), 'MarkerSize', 12, 'LineWidth', 2)
end
plot(tPlot(end) * pi_dt,iterpi_M(end,1), 's', 'MarkerSize', 12, 'LineWidth', 2)
xlabel('$\pi_t$','Interpreter','Latex')
ylabel('$\pi_{\widetilde{M}_{0}}$','Interpreter','Latex')
legend('Combination $1$', 'Combination $2$', 'Combination $3$', 'Combination $4$', 'Combination $5$', 'Non-dimensional', 'location','best', 'FontSize', 14,'Interpreter','Latex')
ylim([0 1])
set(gcf, 'Color', 'w');
set(axes7,'FontSize', 16)
% plot moment pi_M_1 against pi_t for all combinations
figure(8)
axes8 = axes('Parent', figure(8));
hold on
for j = 1:5
    tPlot = itert(:,j);
    iterM1Plot = iterpi_M1(1:stepsize:size(iterpi_M1,1),j);
    plot(tPlot(1:stepsize:size(tPlot,1)) * pi_dt,iterM1Plot, marker(j), 'MarkerSize', 12,  'LineWidth', 2)
end
tPlot = iterpi_t(:);
pi_M1Plot = iterpi_M(1:stepsize:size(iterpi_M,1),2);
plot(tPlot(1:stepsize:size(tPlot,1)) * pi_dt,pi_M1Plot, 's', 'MarkerSize', 12, 'LineWidth', 2)
% reset color order to plot last value
set(axes8, 'ColorOrderIndex', 1);
for j = 1:5
    plot(tPlot(end) * pi_dt, iterpi_M1(end,j), marker(j), 'MarkerSize', 12, 'LineWidth', 2)
end
plot(tPlot(end) * pi_dt,iterpi_M(end,2), 's', 'MarkerSize', 12, 'LineWidth', 2)
xlabel('$\pi_t$','Interpreter','Latex')
ylabel('$\pi_{\widetilde{M}_{1}}$','Interpreter','Latex')
legend('Combination $1$', 'Combination $2$', 'Combination $3$', 'Combination $4$', 'Combination $5$', 'Non-dimensional', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes8,'FontSize', 16)