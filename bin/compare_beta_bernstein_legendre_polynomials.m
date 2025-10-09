%% Test mixture approximations for nonlinear transformations
% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
close all;

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

%% Approximation
% Order
M = 10;

% Interval size
dt = pi;

% Approximation parameters
q = [M, dt];

% Extract inputs
cbeta = [
    0.0001
    0.0017
    0.0220
    0.0509
    0.0842
    0.0005
    0.1829
    0.1604
    0.0001
    0.4801
    0.0173];

cbeta = cbeta(1:M+1);

theta = [];

% Change to Bernstein basis functions
cbernstein = (M+1)/dt*cbeta;

% Basis change matrices (beta to Legendre basis functions)
Mat    = legendre_to_bernstein(q);
Matinv = bernstein_to_legendre(q);
Lambda = beta_to_legendre     (q);

% Change to Legendre basis functions
% clegendre = Mat\cbernstein;
clegendre  = Matinv*cbernstein;
clegendre2 = Lambda*cbeta;

% Define points
t = linspace(0, dt, 1e2);

% Evaluate basis functions
lbeta      = evaluate_beta_basis_functions     (t, theta, q);
lbernstein = evaluate_bernstein_basis_functions(t, theta, q);
llegendre  = evaluate_legendre_basis_functions (t, theta, q);

% Evaluate approximate pdf
alphabeta      = lbeta'     *cbeta;
alphabernstein = lbernstein'*cbernstein;
alphalegendre  = llegendre' *clegendre;
alphalegendre2 = llegendre' *clegendre2;

%% Visualize
% Create figure
figure(1);

% Plot true and approximate pdf
plot(t, alphabeta);      hold on;
plot(t, alphabernstein);
plot(t, alphalegendre); 
plot(t, alphalegendre2); hold off;

% Axis limits
xlim([0, dt]);

% Legend
legend('beta', 'Bernstein', 'Legendre', 'Location', 'NorthWest');