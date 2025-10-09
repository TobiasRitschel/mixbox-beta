function [sn, wn] = chebyshevGaussLobattoNodesWeightsAndWeightFunction(M)
% CHEBYSHEVGAUSSLOBATTONODESWEIGHTSANDWEIGHTFUNCTION Evaluate the Chebyshev
% Gauss-Lobatto nodes, weight, and weight function.
%
% SYNOPSIS:
%   [sn, wn] = chebyshevGaussLobattoNodesWeightsAndWeightFunction(M)
%
% DESCRIPTION:
% Evaluate the Chebyshev Gauss-Lobatto nodes and weights.
% 
% REQUIRED PARAMETERS:
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   sn   - Chebyshev Gauss-Lobatto nodes                    (dimension: M+1)
%   wn   - Chebyshev Gauss-Lobatto weights                  (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussRadauNodesWeightsAndWeightFunction
%           legendreGaussNodesWeightsAndWeightFunction
%           legendrePolynomials
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
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

%% Compute Chebyshev Gauss-Lobatto collocation nodes
% Index
n = (0:M)';

% Collocation nodes
sn = cos(n*pi/M);

% Collocation weights
wn = pi/M*ones(M+1, 1);
wn([1, end]) = pi/(2*M);