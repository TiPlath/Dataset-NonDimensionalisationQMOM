%% Compute kernels for constant growth, aggregation and breakage
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
% 
% Computes the source terms of the population balance equation for
% constant growth, breakage and aggregation kernels. Different kernels can
% be implemented here or in a newly created file.
% 
% INPUT:   xi_alpha         nodes of the N-disperse particle size distribution
%          w_alpha          weights of the N-disperse particle size distribution
%          dotV             growth scale vector for each node
%          a                aggregation scale vector for each node
%          beta             breakage scale for each node
%          b_alpha          fragmentation distribution function for each node
% 
% OUTPUT:  S_V              Vector containing the 2N_delta solutions to the
%                           source terms, i.e. dM_k/dt = S_V

function [S_V] = GrowthAggregationBreakage(xi_alpha,w_alpha,G,a,beta,b_alpha)
    % number of weights
    n = length(w_alpha);
    %% Assemble source term
    % assemble matrices A1, A2 (left hand side of equation)
    A1=zeros(2*n,n);
    A2=zeros(2*n,n);
    % specific for aggregation (right hand side of equation)
    A3=zeros(2*n,n);
    A4=zeros(2*n,n);
    % specific for breakage (right hand side of equation)
    A5=zeros(2*n,n);
    % fragmentation distribution function
    B1=zeros(2*n,n);
    for k=0:2*n-1
        % volume-based growth
        A1(k+1,:)= (1-k)*xi_alpha.^k;
        A2(k+1,:)= k*xi_alpha.^(k-1);
        % volume-based aggregation        
        for i = 1:n
            for r = 1:n
                A3(k+1,i) = A3(k+1,i) + (((xi_alpha(i) + xi_alpha(r))^k) * w_alpha(i) * w_alpha(r));
                A4(k+1,i) = A4(k+1,i) + ((xi_alpha(i)^k) * w_alpha(i) * w_alpha(r));
            end
        end
        % volume-based breakage
        A5(k+1,:) = xi_alpha.^k;
        B1(k+1,:) = b_alpha(xi_alpha,k);   
    end 
    %% compute source terms
    % Aggregation
    Ag = 0.5*(A3*a) - A4*a;
    % Breakage 
    Br = B1*diag(w_alpha)*beta - A5*diag(w_alpha)*beta;
    % Growth
    Gr =  A2*diag(w_alpha)*G;
    % merge breakage, aggregation and growth contribution to source term
    S_V = Ag + Br + Gr;
end