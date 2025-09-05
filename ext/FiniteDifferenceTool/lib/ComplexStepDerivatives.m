function dfFD = ComplexStepDerivatives(fun, x, epsilon, varargin)
% ComplexStepDerivatives  Currently undocumented
%
% REFERENCES
% [1] Al-Mohy, A.H., Higham, N.J., 2010. The complex step approximation to the
% Fréchet derivative of a matrix function. Numerical Algorithms 53, pp. 133-148.
% DOI: 10.1007/s11075-009-9323-y.
% [2] Lai, K.-L., Crassidis, J.L., 2008. Extensions of the first and second
% complex-step derivative approximations. Journal of Computational and Applied
% Mathematics 219(1), pp. 276-293. DOI: 10.1016/j.cam.2007.07.026.
% 
% LINKS
% [A] https://mdolab.engin.umich.edu/wiki/guide-complex-step-derivative-approximation
% [B] https://faculty.math.illinois.edu/~hirani/cbmg/precision.html#:~:text=The%20smallest%20representable%20number%20in,than%201.1102230246316%C3%9710%E2%88%9216.

% Evaluate function
f = feval(fun, x, varargin{:});

% Dimensions
nx = numel(x);
nf = numel(f);

%% First order derivatives
dfFD = zeros(nf, nx);
for j = 1:nx
    % Copy vector
    xp = x;

    % Pertubation
    xp(j) = xp(j) + 1i*epsilon;
    
    % Perturbed function evaluation
    fp = feval(fun, xp, varargin{:});
    
    % Approximation
    dfFD(:, j) = imag(fp(:) / epsilon);
end