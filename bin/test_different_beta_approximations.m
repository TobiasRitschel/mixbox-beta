%% Compare different beta approximations
% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
% close all;

% Remove added paths
restoredefaultpath;

% Restore default settings
reset(groot);

%% Add paths
% Add path to source code
addpath(genpath(fullfile(pwd, '../src')));

%% Formatting
% Save figures?
SaveFigures = true;

% Format of figures
Format = 'eps';

% Font size
fs = 16;

% Line width
lw = 3;

% Set default font size
set(groot, 'DefaultAxesFontSize',   fs);

% Set default line widths
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

% Set default interpreter
set(groot, 'DefaultTextInterpreter',                'Latex');
set(groot, 'DefaultColorbarTickLabelInterpreter',   'Latex');
set(groot, 'DefaultAxesTickLabelInterpreter',       'Latex');
set(groot, 'DefaultLegendInterpreter',              'Latex');

% Set default renderer (otherwise, eps figures can become pixelated in Latex)
set(groot, 'DefaultFigureRenderer', 'Painters');

% Set default units
set(groot, 'DefaultFigureUnits', 'Inches');

% Set default figure size and position
set(groot, 'DefaultFigurePosition', get(groot, 'DefaultFigurePosition'));

%% Initialize
% Maximum delay
dt = 2;

% Midpoint
tm = 1.5;

% True distribution
alpha = @(t) 2*t/(tm*dt).*(t <= tm) + 2*(dt - t)./(dt.*(dt - tm)).*(t > tm);

%% Estimate coefficients
% Number of points
N = 1e2;

% Order
M = 20;

% QP estimate
[clsb,        lambdalsb ] = least_squares_estimate_qp           (alpha, dt, M, N, @least_squares_qp_beta);
[clsl,        lambdalsl ] = least_squares_estimate_qp           (alpha, dt, M, N, @least_squares_qp_legendre);
[clsl2,       lambdalsl2] = least_squares_estimate_qp_legendre  (alpha, dt, M);
[cint, muint, lambdaint ] = least_squares_interpolation_estimate(alpha, dt, M   );

% Times
t = linspace(0, dt, 1e4);

% Approximation parameters
q     = [M, dt];
theta = [];

% Evaluate basis functions
lbeta     = evaluate_beta_basis_functions    (t, theta, q);
llegendre = evaluate_legendre_basis_functions(t, theta, q);

% Conversion matrix
A = beta_to_legendre(q);

% Evaluate approximate kernel
alphahatlsb  = lbeta'    *clsb;
alphahatlsl  = llegendre'*clsl;
alphahatlsl2 = llegendre'*clsl2;
alphahatint  = lbeta'    *cint;

%% Visualize
% Create figure
figure(1); clf;

% Select subplot
subplot(211);

% Plot true and approximate pdf
plot(t, alpha(t)); hold on;
plot(t, alphahatlsb);
plot(t, alphahatlsl);
plot(t, alphahatlsl2);
plot(t, alphahatint); hold off;

% Axis limits
xlim([0, dt]);

% Legend
legend('True', 'Least-squares (beta)', 'Least-squares (Legendre)', 'Least-squares (Legendre)', 'Interpolation', 'Location', 'NorthWest');

% Select subplot
subplot(212);

% Color index
set(gca, 'ColorOrderIndex', 2);

% Keep adding
hold on;

% Plot true and approximate pdf
semilogy(t, abs(alpha(t)' - alphahatlsb ));
semilogy(t, abs(alpha(t)' - alphahatlsl ));
semilogy(t, abs(alpha(t)' - alphahatlsl2));
semilogy(t, abs(alpha(t)' - alphahatint )); hold off;

% Axis limits
xlim([0, dt]);

%% Functions
function [c, lambda] = least_squares_estimate_qp_legendre(alpha, dt, M)
% Legendre nodes and weights
[tau, v] = legendreGaussLobattoNodesWeightsAndWeightFunction(M);

% Time points
t = dt/2*(tau + 1)';

% Weights
wbar = dt/2*v;

% Approximation parameters
q = [M, dt];
theta = [];

% Conversion matrix
A = beta_to_legendre(q);

% Evaluate basis functions
L      = evaluate_legendre_basis_functions(t, theta, q);
alphak = alpha(t);

% Quadratic coefficients
H = diag(dt./(2*(0:M) + 1));

% Coefficients
e = A'\ones(M+1, 1);
g = -L*(alphak'.*wbar);

% System matrix and right-hand side
A = [H, e; e', 0];
b = [-g; 1];

% Estimate coefficients
x = A\b;

% Coefficients and Lagrange multiplier
c      = x(1:M+1);
lambda = x(M+2);
end

function [c, lambda] = least_squares_estimate_qp(alpha, dt, M, N, evaluate_least_squares_qp)
% Compute matrix and vector
[H, g, e] = evaluate_least_squares_qp(alpha, dt, M, N);

% System matrix and right-hand side
A = [H, e; e', 0];
b = [-g; 1];

% Estimate coefficients
x = A\b;

% Coefficients and Lagrange multiplier
c      = x(1:M+1);
lambda = x(M+2);
end

function [H, g, ebar] = least_squares_qp_legendre(alpha, dt, M, N)
% Approximation parameters
q = [M, dt];

% Conversion matrix
A = beta_to_legendre(q);

% Evaluate coefficients based on beta basis functions
[H, g, e] = least_squares_qp(alpha, dt, M, N, @evaluate_legendre_basis_functions);

% Coefficients for Legendre basis
ebar = A'\e;
end

function [H, g, e] = least_squares_qp_beta(alpha, dt, M, N)
    [H, g, e] = least_squares_qp(alpha, dt, M, N, @evaluate_beta_basis_functions);
end

function [H, g, e] = least_squares_qp(alpha, dt, M, N, evaluate_basis_functions)
% Time points
t = dt/2*(legendreGaussLobattoNodesWeightsAndWeightFunction(N) + 1)';

% Approximation parameters
q     = [M, dt];
theta = [];

% Evaluate basis functions
lk     = evaluate_basis_functions(t, theta, q);
alphak = alpha(t);

% Linear coefficients
g = -sum(alphak.*lk, 2);
e = ones(M+1, 1);

% Initialize
H = 0;

for k = 1:N+1
    % Add to matrix
    H = H + lk(:, k)*lk(:, k)';
end
end

function [c, mu, lambda] = least_squares_interpolation_estimate(alpha, dt, M)
% Time points
t = dt/2*(legendreGaussLobattoNodesWeightsAndWeightFunction(M) + 1)';

% Approximation parameters
q     = [M, dt];
theta = [];

% Evaluate basis functions
lk     = evaluate_beta_basis_functions(t, theta, q);
alphak = alpha(t);

% Matrix
A = [eye(M+1),    lk, ones(M+1, 1);
    lk',          zeros(M+1, M+2);
    ones(1, M+1), zeros(  1, M+2)];

% Right-hand side
b = [zeros(M+1, 1); alphak'; 1];

% Estimate coefficients
x = A\b;

% Coefficients and Lagrange multiplier
c      = x(          1:M+1);
mu     = x(   M+1 + (1:M+1));
lambda = x(2*(M+1) + 1     );
end

function [c, mu, lambda] = least_squares_interpolation_estimate_explicit_DEPRECATED(alpha, dt, M)
% Time points
t = dt/2*(legendreGaussLobattoNodesWeightsAndWeightFunction(M) + 1)';

% Approximation parameters
q     = [M, dt];
theta = [];

% Evaluate basis functions
lk     = evaluate_beta_basis_functions(t, theta, q);
alphak = alpha(t);

% Matrices and vectors
I = eye(M+1);
e = ones(M+1, 1);

A = lk;
b = alphak';

% Auxiliary quantity
aux1 = (A'*A)\b;
aux2 = (A'*A)\A';

% Lagrange multipliers
lambda = (1 - e'*A*aux1)/(e'*(A*aux2 - I)*e);
mu     = -aux1 - aux2'*e*lambda;

% Coefficients
c = -A*mu - e*lambda;
end