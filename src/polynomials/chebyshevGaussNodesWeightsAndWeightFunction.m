function [sn, wn] = chebyshevGaussNodesWeightsAndWeightFunction(M)
% CHEBYSHEVGAUSSNODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Chebyshev
% Gauss nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = chebyshevGaussNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Chebyshev Gauss nodes and weights.
% 
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Chebyshev Gauss nodes                    (dimension: M+1)
%   wn   - Chebyshev Gauss weights                  (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussRadauNodesWeightsAndWeightFunction
%           legendreGaussNodesWeightsAndWeightFunction
%           legendrePolynomials
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
%           chebyshevGaussLobattoNodesWeightsAndWeightFunction
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

%% Compute Chebyshev Gauss collocation nodes
% Index
n = (0:M)';

% Collocation nodes
sn = cos((2*n+1)*pi/(2*M+2));

% Collocation weights
wn = pi/(M+1)*ones(M+1, 1);