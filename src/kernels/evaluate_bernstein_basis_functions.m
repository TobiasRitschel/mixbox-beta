function l = evaluate_bernstein_basis_functions(t, ~, q)

% Extract parameters
M  = q(1);
dt = q(2);

% Parametrize time
s = t/dt;

% Orders
m = (0:M)';

% Coefficients
b = arrayfun(@nchoosek, M*ones(M+1, 1), m);

% Evaluate basis function
l = b.*s.^m.*(1 - s).^(M - m);