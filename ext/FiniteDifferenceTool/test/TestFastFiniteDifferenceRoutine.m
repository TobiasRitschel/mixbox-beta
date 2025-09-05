%% Test finite difference routine
clc; clear all; close all

%% Load library and utility routines
% Load library
run('../LoadLibrary');

%% Initialization
fun = @FiniteDifferenceTestFunction;

% Test point
x = [1; 2];

% Finite differences
NoTests = 1e2;
epsilon = 10.^(linspace(-8, 2, NoTests));

% Error function
ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:)) ./ (1 + abs(dfFD(:)));

% Carry out finite difference test
FastFiniteDifferenceTest(fun, x, epsilon, true, ErrorFunction);