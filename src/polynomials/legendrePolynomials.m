function [LM, dLM] = legendrePolynomials(z, M)
% LEGENDREPOLYNOMIALS Evaluate Legendre polynomials and their derivatives.
%
% SYNOPSIS:
%   [L, dL] = legendrePolynomials(z, M)
%
% DESCRIPTION:
% Evaluate Legendre polynomials and their first order derivatives.
% 
% REQUIRED PARAMETERS:
%   z - spatial coordinate
%   M - order of polynomial
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   l   - Legendre polynomials                              (dimension: M+1)
%   dl  - First order derivatives of Legendre polynomials   (dimension: M+1)
%
% DEPENDENCIES:
%
% See also  legendreGaussLobattoNodesWeightsAndWeightFunction
%           legendreGaussNodesWeightsAndWeightFunction
%           chebyshevGaussLobattoNodesWeightsAndWeightFunction
%           chebyshevGaussNodesWeightsAndWeightFunction
%           chebyshevGaussRadauNodesWeightsAndWeightFunction
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
%  info@diamatica.com
%  tobk@dtu.dk
% athre@dtu.dk
%  jbjo@dtu.dk
% 
% AUTHORS
% Tobias K. S. Ritschel
% Asbjørn Thode Reenberg
% John Bagterp Jørgensen

%% Evaluate Legendre polynomial
if(M == 0)
    LM  = 1;
    dLM = 0;
elseif(M == 1)
    LM  = z;
    dLM = 1;
else
    % Initialize polynomials
    LMm2 = 1;
    LMm1 = z;
    
    % Initialize derivatives of polynomials
    dLMm2 = 0;
    dLMm1 = 1;
    
    for k = 2:M
        % Compute polynomial and its derivative
        LM  = (2*k - 1)/k*z.*LMm1 - (k - 1)/k*LMm2;
        dLM = dLMm2 + (2*k - 1)*LMm1;
        
        % Update "previous" polynomials and their derivatives
        LMm2 = LMm1;
        LMm1 = LM;
        
        dLMm2 = dLMm1;
        dLMm1 = dLM;
    end
end