function [alpha, dalpha, l] = evaluateKernel(t, c, dt, varargin)

% Compute gradient and Hessian?
ComputeGradient = (nargout > 1);

% Type of basis function
type = 'beta';
if(nargin > 3), type = lower(varargin{1}); end

% Stable computations?
stable = 'yes';
if(nargin > 4), stable = varargin{2}; end

% Maximal order
M  = numel(c)-1;

switch type
    case 'beta'
        % Evaluate beta basis functions
        l = evaluate_beta_basis_functions(t, [], [M, dt], stable);

    case 'legendre'
        % Evaluate beta basis functions
        l = evaluate_legendre_basis_functions(t, [], [M, dt], stable);

    case 'bernstein'
        % Evaluate beta basis functions
        l = evaluate_bernstein_basis_functions(t, [], [M, dt], stable);
end

% Mixture kernel
alpha = c*l;

if(ComputeGradient)
    % Allocate memory
    dalpha = zeros(1, numel(t), M+1);

    for idx = 1:M+1
        % Derivative of kernel wrt. coefficient
        dalpha(:, :, idx) = l(idx, :);
    end
end