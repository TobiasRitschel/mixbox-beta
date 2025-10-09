function l = evaluate_beta_basis_functions(x, ~, q, varargin)
% Extract parameters
M  = q(1);
dt = q(2);

% Stable computations?
stable = 'yes';
if(nargin > 3), stable = varargin{1}; end

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
l = b.*x.^m.*(dt - x).^(M - m);