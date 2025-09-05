function [c, obj, exitflag, output] = identify_kernel(t, M, alpha, dt, varargin)

% Options
opts = [];
if(nargin > 4), opts = varargin{1}; end

% fmincon options
fmincon_opts = optimoptions('fmincon');
if(~isempty(opts) && isfield(opts, 'fmincon_opts') && ~isempty(opts.fmincon_opts))
    fmincon_opts = opts.fmincon_opts;
end

% Matlab routines with double kernel computations
objective_gradient_routine = @least_squares_kernel_objective;
hessian_routine            = @least_squares_kernel_objective_hessian;

% Evaluate true kernel
alphameas = alpha(t, varargin{2:end});

% Overwrite specific fmincon settings
fmincon_opts = optimoptions(fmincon_opts, ...
    'SpecifyObjectiveGradient', true,          ...
    'HessianFcn',               hessian_routine);

%% Estimate kernel parameters
% Initial guess of interpolation coefficients and forgetting rate
c0 = ones(1, M+1)/(M+1);

% Inequality constraint system matrix and right-hand side
A = [];
B = [];

% Equality constraint system matrix and right-hand side
Aeq = ones(1, M+1);
Beq = 1;

% Upper and lower bounds
ub =  ones(M+1, 1);
lb = zeros(M+1, 1);

% Nonlinear constraint function
nonlcon = [];

% Estimate the interpolation coefficients
[c, obj, exitflag, output, lambda, grad, hessian] = ...
    fmincon(objective_gradient_routine, ...
    c0, A, B, Aeq, Beq, lb, ub, nonlcon, fmincon_opts, ...
    t, alphameas, dt, opts); %#ok