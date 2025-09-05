function hess = least_squares_kernel_objective_hessian(c, ~, t, ~, dt, opts)

% Number of parameters
ntheta = numel(c);

% Reshape
c = c(:)';

% Scaling factors
if(~isempty(opts) && isfield(opts, 'rho')), rho = opts.rho; else, rho = 1; end

% Indices of points used in quadrature
idx = 1:numel(t)-1;

% Choose the left points in the integral discretization
tl = t(idx);

% Compute the time intervals
dti = diff(t);

% Evaluate approximate kernel
[~, dalphahat] = evaluateKernel(tl, c, dt);

% Derivative of deviation
de = -dalphahat;

% Allocate memory
hess = zeros(ntheta, ntheta);

for i = 1:ntheta
    for j = 1:ntheta
        % Hessian element
        hess(i, j) = rho*(de(:, :, j).*dti)*de(:, :, i)';
    end
end