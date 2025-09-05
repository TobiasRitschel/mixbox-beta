function [alpha, dalpha, l] = evaluateKernel(t, c, dt)

% Compute gradient and Hessian?
ComputeGradient = (nargout > 1);

%% Mixed Beta kernel
% Maximal order
M  = numel(c)-1;

% Orders
m = (0:M)';

% Coefficients
b = (M+1)/dt^(M+1)*arrayfun(@nchoosek, M*ones(M+1, 1), m);

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