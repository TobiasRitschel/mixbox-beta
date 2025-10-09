function [sn, wn] = legendreGaussNodesWeightsAndWeightFunction(M)
% LEGENDREGAUSSNODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Legendre
% Gauss nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = legendreGaussNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Legendre Gauss nodes and weights.
%
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Legendre Gauss nodes                         (dimension: M+1)
%   wn   - Legendre Gauss weights                       (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  chebyshevGaussLobattoNodesWeightsAndWeightFunction
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
%           chebyshevGaussNodesWeightsAndWeightFunction
%           legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussRadauNodesWeightsAndWeightFunction
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

%% Compute Legendre Gauss collocation nodes
% Options
opts = optimoptions('fsolve',   ...
    'Display',      'off',      ...
    'TolFun',       eps/2);

if(M == 0)
    % End points
    sn(1) = 0;
    wn(1) = 0;
elseif(M == 2)
    % End points
    sn(1) = -sqrt(1/3);
    wn(1) =  1;
    sn(2) = -sn(1);
    wn(2) =  wn (1);
else
    % Allocate memory
    sn = zeros(M+1, 1);
    wn = zeros(M+1, 1);
    
    for n = 0:(M + 1)/2 - 1
        % Initial guess (Chebyshev Gauss collocation nodes)
        sn0 = -cos((2*n + 1)/(2*M + 2)*pi);
        
        % Residual function
        Res = @(s) legendrePolynomials(s, M+1);
        
        % Solve for the root
        sn(n+1) = fsolve(Res, sn0, opts);
        
        % Copy root
        sn(M+1-n) = -sn(n+1);
        
        % Evaluate Legendre polynomial
        [~, dLMp1] = legendrePolynomials(sn(n+1), M+1);
        
        % Weights
        wn(n+1) = 2/((1 - sn(n+1)^2)*dLMp1^2);
        
        % Copy weight
        wn(M+1-n) = wn(n+1);
    end
    
    if(mod(M, 2) == 0)
        % Middle point
        sn(M/2+1) = 0;
        
        % Evaluate Legendre polynomial
        [~, dLMp1] = legendrePolynomials(sn(M/2+1), M+1);
        
        % Middle weight
        wn(M/2+1) = 2/dLMp1^2;
    end
end