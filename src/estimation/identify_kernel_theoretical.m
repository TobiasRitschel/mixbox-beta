function c = identify_kernel_theoretical(alphabeta, dt, M, varargin)

% Load option
option = 'integral';
if(nargin > 3), option = varargin{1}; end

switch lower(option)
    case 'integral'
        % Integral
        beta = alphabeta;

        % Time step
        dti = dt/(M+1);

        % Create time grid (used to evaluate coefficients)
        tbar = 0:dti:dt;

        % Evaluate coefficients
        c = diff(beta(tbar));
    case 'function'
        % Function
        alpha = alphabeta;

        % Create time grid (used to evaluate coefficients)
        tbar = linspace(0, dt, M+1);

        % Evaluate coefficients
        c = dt*alpha(tbar)/(M+1);
end