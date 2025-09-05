% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
close all;

% Remove added paths
restoredefaultpath;

% Reset default settings
reset(groot);

%% Add paths
% Add path to source code
run('../load_library.m');

%% Formatting
% Save figures?
SaveFigures = true;

% Format of figures
Format = 'eps';

% Font size
fs = 26;

% Line width
lw = 4;

% Marker size
ms = 20;

% Face alpha
fa = 0.9;

% Get default figure position
defaultFigurePosition = get(groot, 'DefaultFigurePosition');

% Set default font size
set(groot, 'DefaultAxesFontSize',   fs);
set(groot, 'DefaultTextFontSize',   fs);

% Set default line widths
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

% Set default marker size
set(groot, 'DefaultLineMarkerSize', ms);

% Set default interpreter
set(groot, 'DefaultTextInterpreter',            'Latex');
set(groot, 'DefaultAxesTickLabelInterpreter',   'Latex');
set(groot, 'DefaultLegendInterpreter',          'Latex');

% Increase default figure position
set(groot, 'DefaultFigurePosition', defaultFigurePosition.*[1, 1, 2*0.8, 2*0.6]);

%% Kernel
% Domain
dt = pi;

% Radius
R = dt/2;

% Kernel
alpha = @(t) 2*sqrt(R^2 - (t - R).^2)/(pi*R^2);

%% Approximate kernel
% Order
M = 2e1;

% Number of measurements
N = 200;

% Measurement points
tmeas = linspace(0, dt, N+1);

% Identify kernel for fixed order
[chat, obj, exitflag, output] = identify_kernel(tmeas, M, alpha, dt);

% Report results
fprintf('\n');
fprintf('c%-3d = %7.3f\n', [0:M; chat]);

%% Visualize kernel
% Times
tkern = linspace(0, dt, 1e3);

% Estimated kernel
[alphahat, ~, alphahat_subkernels] = evaluateKernel(tkern, chat, dt);

% Axis tick labels for bar chart of estimated variables
XTickLabels = cell(1, M+1);

for l = 1:M+1
    XTickLabels{l} = "$c_{" + num2str(l-1) + "}$";
end

%% True kernel
% Select subplot
subplot(221);

% Visualize kernel
plot(tkern, alpha(tkern), 'DisplayName', 'True', 'Color', [0, 0.4470, 0.7410]);

% Axis limits
xlim(tkern([1, end]));

% Title
title('True density');

% Axis handle
h = gca;

% Change font size
h.FontSize = fs;

%% Approximate kernel
% Select subplot
subplot(223);

% Visualize kernel
plot(tkern, alphahat, 'DisplayName', ['Estimated (M = ', num2str(M), ')']);

% Axis limits
xlim(tkern([1, end]));

% Title
title('Approximate density');

% Axis handle
h = gca;

% Change font size
h.FontSize = fs;

%% Subkernels
% Select subplot
subplot(222);

% Visualize kernel
plot(tkern, alphahat_subkernels, 'DisplayName', ['Estimated (M = ', num2str(M), ')']);

% Axis limits
xlim(tkern([1, end]));

% Title
title('Basis functions');

% Axis handle
h = gca;

% Change font size
h.FontSize = fs;

%% Coefficients
% Select figure
figure(1);

% Select subplot
subplot(224);

% Dummy plot to keep box
bar(NaN);

% Set color index
set(gca, 'ColorOrderIndex', 2);

% Add more plots
hold on;

% Visualize true parameter estimates
b = histogram('BinEdges', (0:M+1)+0.5, 'BinCounts', chat, 'FaceAlpha', fa);

% Stop adding plots
hold off;

% Title
title('Coefficients');

% Axis limits
xlim([-0.75, M+2.75]);

% Axis ticks
xticks(1:10:M+1);

% Axis tick labels
xticklabels(XTickLabels(1:10:M+1));

% Remove bar edges
for i = 1:numel(b), b(i).EdgeColor = 'none'; end

% Remove tick lines
h = gca;
h.XAxis.TickLength = [0, 0];

% Axis handle
h = gca;

% Change font size
h.FontSize = fs;

% Save plot
SavePlot('erlang_mixture_approximation_overview', SaveFigures, Format);

%% Error
% Increase default figure position
set(groot, 'DefaultFigurePosition', defaultFigurePosition.*[1, 1, 0.8, 0.6]);

% Create figure
figure(3);

% Visualize kernel
semilogy(tkern, abs(alpha(tkern) - alphahat), 'DisplayName', 'True', 'Color', [0, 0.4470, 0.7410]);

% Axis limits
xlim(tkern([1, end]));

% Title
title('Approximation error');

% Axis handle
h = gca;

% Change font size
h.FontSize = fs;

% Save plot
SavePlot('erlang_mixture_approximation_error', SaveFigures, Format);