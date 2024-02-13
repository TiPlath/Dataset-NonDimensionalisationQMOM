%% Quadrature Method of Moments non-dimensionalized
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% This script utilizes the quadrature method of moments for combined
% aggregation, breakage and growth mechanisms. The kernels can be adjusted if needed
% It reads in a PSD and evolves its set of moments in time. One can compute 
% physical values from these moments to get information about the evolved 
% particle size distribution.
% 
% INPUT: A .csv file for the function ReadData, containing information of 
%        a particle size distribution, i.e. diameters (column 1) and
%        probabilities/densities (column 2).
% 
% OUTPUT: iterM0    A vector containing the zeroth moment for all
%                   combinations and timesteps computed.
%         iterM1    A vector containing the first moment for all
%                   combinations and timesteps computed.

%% MATLAB options
clear all
close all
clc
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% read custom Distribution
combinations = 5;
% timestep
dt=0.1;
% max time
tmax=10;
% Number of time steps
nt=round(tmax/dt)+1;
% store M0,M1,t for every iteration with different aggregation rates
iterM0 = zeros(nt,combinations);
iterM1 = zeros(nt,combinations);
itert = zeros(nt,combinations);
for j = 1:combinations
    % call function to read distribution with units. In our case we read in a
    % volume-based number distribution (n_V).
    n_V = ReadData('Data/InitialLactoseMCCPVP-n_V.csv',1,0,1);
    %% Define variables (pre-processing)
    % Number of dirac-delta distributed classes (weights and nodes, 1,...,25)
    N_delta=1;
    % Maximum number of moments
    mMax = 2*N_delta;
    % Time
    t=(0:nt);
    % A matrix storing the momenta of every time step
    iterM=zeros(nt,2*N_delta);
    % store weights and nodes of every iteration
    iterV_alpha = zeros(nt,N_delta);
    iterW_alpha = zeros(nt,N_delta);    
    %%  Kernel values
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
    for l = 1:size(beta)
        beta(l,:) = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
    end
    % symmetric fragmentation for n_V
    b_alpha = @(V_alpha,k) 2^((1-k))*V_alpha.^k;
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
%% Post-Processing
% define stepsize for plotting to have less markers
stepsize = nt*dt;
% plot moment M_0 against t for all combinations
marker = ['x','o','+','d','v'];
figure(7)
axes7 = axes('Parent', figure(7));
hold on
for j = 1:5
    tPlot = itert(:,j);    
    iterM0Plot = iterM0(1:stepsize:size(iterM0,1),j)/iterM0(1,j);
    plot(tPlot(1:stepsize:size(tPlot,1)) * dt,iterM0Plot, marker(j), 'MarkerSize', 12,  'LineWidth', 2)
end
% reset color order to plot last value
set(axes7, 'ColorOrderIndex', 1);
for j = 1:5
    plot(tPlot(end) * dt, iterM0(end,j)/iterM0(1,j), marker(j), 'MarkerSize', 12, 'LineWidth', 2)
end
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$\hat{M}_{0}$ [-]','Interpreter','Latex')
legend('Combination $1$', 'Combination $2$', 'Combination $3$', 'Combination $4$', 'Combination $5$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
ylim([0 1])
set(gcf, 'Color', 'w');
set(axes7,'FontSize', 16)
% plot moment M_1 against t for all combinations
figure(8)
axes8 = axes('Parent', figure(8));
hold on
for j = 1:5
    tPlot = itert(:,j);
    iterM1Plot = iterM1(1:stepsize:size(iterM1,1),j)/iterM1(1,j);
    plot(tPlot(1:stepsize:size(tPlot,1)) * dt, iterM1Plot, marker(j), 'MarkerSize', 12,  'LineWidth', 2)
end
% reset color order to plot last value
set(axes8, 'ColorOrderIndex', 1);
for j = 1:5
    plot(tPlot(end) * dt, iterM1(end,j)/iterM1(1,j), marker(j), 'MarkerSize', 12, 'LineWidth', 2)
end

xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$\hat{M}_{1}$ [-]','Interpreter','Latex')
legend('Combination $1$', 'Combination $2$', 'Combination $3$', 'Combination $4$', 'Combination $5$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes8,'FontSize', 16)