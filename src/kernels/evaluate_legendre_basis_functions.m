function l = evaluate_legendre_basis_functions(t, ~, q)

% Extract parameters
M  = q(1);
dt = q(2);

% Parametrize time
s = t/dt;

% Allocate memory
l = zeros(M+1, numel(t));

% Initial polynomials
l(1, :) =        1;
l(2, :) = (2*s - 1);

for m = 2:M
    % Evaluate Legendre polynomials recursively
    l(m+1, :) = ((2*m - 1)*(2*s - 1).*l(m, :) - (m - 1)*l(m-1, :))/m;
end