%% Compute moments of a distribution
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Computes the initial moments from a given function and returns them in a
% Vector M. The computed moments are raw moments according to the formula
% M_k = \int_x (x^r * f(x))
% 
% INPUT:   f             probability/density values of the PSD
%          x             Domain of the PSD
%          mMax          Maximum number of moments 
% 
% OUTPUT:  M             Matrix  with desired Moments of the current distribution
%                        with the range from 2...50 moments (50x25 Matrix at max)

function M = ComputeMoments(x, f, mMax)
% Initialize moment matrix
M = zeros(mMax,mMax);
% Zeroth order raw Moment
M(1,1:mMax) = sum(f(1:end-1).*diff(x));
% Higher order  Moments
for N = 2:mMax
    for r = 1:N-1
        M(r+1,N) = sum(f(1:end-1).*diff(x).*(x(1:end-1).^r));
    end
end
end