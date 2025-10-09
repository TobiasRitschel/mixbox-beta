function [sn, wn] = legendreGaussLobattoNodesWeightsAndWeightFunction(M)
% LEGENDREGAUSSLOBATTONODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Legendre
% Gauss-Lobatto nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = legendreGaussLobattoNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Legendre Gauss-Lobatto nodes and weights.
%
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Legendre Gauss-Lobatto nodes                    (dimension: M+1)
%   wn   - Legendre Gauss-Lobatto weights                  (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  chebyshevGaussLobattoNodesWeightsAndWeightFunction
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
%           chebyshevGaussNodesWeightsAndWeightFunction
%           legendreGaussRadauNodesWeightsAndWeightFunction
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

%% Compute Legendre Gauss-Lobatto collocation nodes
% Options
opts = optimoptions('fsolve',   ...
    'Display',      'off',      ...
    'TolFun',       1e-15);

% Allocate memory
sn = zeros(M+1, 1);
wn = zeros(M+1, 1);

if(M == 1)
    % End points
    sn(1) = -1;
    wn(1) =  1;
    sn(2) =  1;
    wn(2) =  1;
else
    % End points
    sn(  1) = -1;
    wn(  1) =  2/(M*(M+1));
    sn(M+1) =  1;
    wn(M+1) = wn(1);
    
    for n = 1:(M + 1)/2 - 1
        % Initial guess
        sn0 = -cos((n + 0.25)*pi/M - 3/(8*M*pi*(n + 0.25)));
        
        % Residual function
        Res = @(s) legendrePolynomials(s, M+1) - legendrePolynomials(s, M-1);
        
        % Solve for the root
        sn(n+1) = fsolve(Res, sn0, opts);
        
        % Copy root
        sn(M+1-n) = -sn(n+1);
        
        % Evaluate Legendre polynomial
        LM = legendrePolynomials(sn(n+1), M);
        
        % Weights
        wn(n+1) = 2/(M*(M+1)*LM^2);
        
        % Copy weight
        wn(M+1-n) = wn(n+1);
    end
    
    if(mod(M, 2) == 0)
        % Middle point
        sn(M/2+1) = 0;
        
        % Evaluate Legendre polynomial
        LM = legendrePolynomials(sn(M/2+1), M);
        
        % Middle weight
        wn(M/2+1) = 2/(M*(M+1)*LM^2);
    end
end