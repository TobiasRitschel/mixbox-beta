function [sn, wn] = chebyshevGaussRadauNodesWeightsAndWeightFunction(M)
% CHEBYSHEVGAUSSRADAUNODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Chebyshev
% Gauss-Radau nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = chebyshevGaussRadauNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Chebyshev Gauss-Radau nodes and weights.
% 
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Chebyshev Gauss-Radau nodes                    (dimension: M+1)
%   wn   - Chebyshev Gauss-Radau weights                  (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussRadauNodesWeightsAndWeightFunction
%           legendreGaussNodesWeightsAndWeightFunction
%           legendrePolynomials
%           chebyshevGaussLobattoNodesWeightsAndWeightFunction
%           chebyshevGaussNodesWeightsAndWeightFunction
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

%% Compute Chebyshev Gauss-Radau collocation nodes
% Index
n = (0:M)';

% Collocation nodes
sn = cos(2*n*pi/(2*M+1));

% Collocation weights
wn = 2*pi/(2*M+2)*ones(M+1, 1);
wn(1) = pi/(2*M+1);