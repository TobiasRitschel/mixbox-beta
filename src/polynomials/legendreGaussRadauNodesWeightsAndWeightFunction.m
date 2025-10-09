function [sn, wn] = legendreGaussRadauNodesWeightsAndWeightFunction(M)
% LEGENDREGAUSSRADAUNODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Legendre
% Gauss-Radau nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = legendreGaussRadauNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Legendre Gauss-Radau nodes and weights.
%
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Legendre Gauss-Radau nodes                    (dimension: M+1)
%   wn   - Legendre Gauss-Radau weights                  (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  chebyshevGaussLobattoNodesWeightsAndWeightFunction
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
%           chebyshevGaussNodesWeightsAndWeightFunction
%           legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussNodesWeightsAndWeightFunction
%           legendrePolynomials
%           lagrangeInterpolation
%           lagrangePolynomials
%           lagrangePolynomialDerivatives
%
% REFERENCES
% [1] D. A. Kopriva, Implementing spectral methods for partial differential
% equations: Algorithms for scientists and engineers. Scientific Computation,
% Springer, 2009.
% [2] C. Canuto, M. Y. Hussaini, A. Quarteroni, and T. A. Zang, Spectral
% methods: Fundamentals in single domains, ser. Scientific Computation.
% Springer, 2006.
%
% CONTACT INFORMATION
% tobk@dtu.dk
%
% AUTHORS
% Tobias K. S. Ritschel

%% Compute Legendre Gauss-Radau collocation nodes
% Options
opts = optimoptions('fsolve',   ...
    'Display',      'none',      ...
    'TolFun',       1e-15);

% Allocate memory
sn = zeros(M+1, 1);
wn = zeros(M+1, 1);

if(M == 0)
    % End points
    sn(1) = 1;
    wn(1) = 2;
else
    for n = 0:M
        % Initial guess (Chebyshev Gauss-Radau collocation nodes)
        sn0 = -cos(2*n/(2*M + 1)*pi);
        
        % Residual function
        Res = @(s) legendrePolynomials(s, M+1) + legendrePolynomials(s, M);
        
        % Solve for the root
        sn(n+1) = fsolve(Res, sn0, opts);
        
        % Evaluate Legendre polynomial
        LM = legendrePolynomials(sn(n+1), M);
        
        % Weights
        switch n
            case 0
                wn(n+1) = 2/(M+1)^2;
            otherwise
                wn(n+1) = (1 - sn(n+1))/((M+2)^2*LM^2);
        end
    end
end
end