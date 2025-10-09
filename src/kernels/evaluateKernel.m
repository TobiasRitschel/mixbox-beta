function [alpha, dalpha, l] = evaluateKernel(t, c, dt, varargin)

% Compute gradient and Hessian?
ComputeGradient = (nargout > 1);

% Stable computations?
stable = 'yes';
if(nargin > 3), stable = varargin{1}; end

%% Mixed Beta kernel
% Maximal order
M  = numel(c)-1;

if(stable)
    % Allocate memory
    b = zeros(M+1, 1);

    % Normalization constants (more stable than nchoosek for large M)
    b(1) = (M+1)/dt^(M+1);
    for m = 1:M
        b(m+1) = b(m)*(M - (m-1))/m;
    end
else
    % Coefficients
    b = (M+1)/dt^(M+1)*arrayfun(@nchoosek, M*ones(M+1, 1), 0:M);
end

% Orders
m = (0:M)';

% Auxiliary variable
l = b.*t.^m.*(dt - t).^(M - m);

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