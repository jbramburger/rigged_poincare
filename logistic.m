% -------------------------------------------------------------------------
% EDMD, mpEDMD, AND riggedDMD FOR THE LOGISTIC MAP
% -------------------------------------------------------------------------
%
% This script reproduces all six panels in the two logistic-map figures:
%
%   logistic_eigenvalues_1.pdf          d = 11
%   logistic_eigenvalues_2.pdf          d = 51
%   logistic_eigenfunction.pdf          d = 51, mpEDMD eigenvalue nearest -1
%   logistic_measure_1.pdf              d = 51, epsilon = 0.01
%   logistic_measure_2.pdf              d = 51, epsilon = 0.75
%   logistic_rigged_eigenfunction.pdf   d = 51, epsilon = 0.75, theta = pi
%
% The trajectory contains 20001 retained states. A delay dictionary of
% dimension d therefore gives 20001-d paired rows: 19990 rows for d = 11
% and 19950 rows for d = 51. Each computation uses uniform empirical
% weights normalized over the rows in its own delay dictionary.
%
% The riggedDMD implementation used here is available from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% Parameters

r = 3.7;
x0 = 0.25;
nTransient = 100;
nRetainedTransitions = 20000;

delayDimensions = [11, 51];
riggedDimension = 51;

riggedEpsilons = [0.01, 0.75];
thetaTarget = pi;
thetaGrid = -pi:0.01:pi;
kernelOrder = 2;

kneadingLength = 40;
criticalPoint = 0.5;

% The PDFs are written to the current directory. Change this value if the
% manuscript expects a separate figures directory.
outputDirectory = '.';

if ~isfolder(outputDirectory)
    mkdir(outputDirectory)
end

%% Generate 20001 retained states

nTotalTransitions = nTransient + nRetainedTransitions;
xAll = zeros(nTotalTransitions + 1, 1);
xAll(1) = x0;

for n = 1:nTotalTransitions
    xAll(n + 1) = logistic_map(xAll(n), r);
end

% Keeping the state at the end of the transient leaves 20001 states and
% hence 20000 one-step pairs.
x = xAll(nTransient + 1:end);
nRetainedStates = numel(x);

assert(nRetainedStates == nRetainedTransitions + 1, ...
    'The retained trajectory must contain 20001 states.');

fprintf('Retained states: %d\n', nRetainedStates);
fprintf('Retained one-step pairs: %d\n', nRetainedStates - 1);
fprintf('Sampled attractor range: [%.16g, %.16g]\n', min(x), max(x));

%% Compute and report the kneading sequence

criticalOrbit = zeros(kneadingLength + 1, 1);
criticalOrbit(1) = criticalPoint;

for n = 1:kneadingLength
    criticalOrbit(n + 1) = logistic_map(criticalOrbit(n), r);
end

kneadingSequence = repmat(' ', kneadingLength, 1);
for n = 1:kneadingLength
    if abs(criticalOrbit(n + 1) - criticalPoint) < 1e-12
        kneadingSequence(n) = 'C';
    elseif criticalOrbit(n + 1) < criticalPoint
        kneadingSequence(n) = 'L';
    else
        kneadingSequence(n) = 'R';
    end
end

fprintf('Kneading sequence: %s\n', kneadingSequence.');

%% Center and normalize the scalar observable

observable = x - mean(x);
observable = observable / std(observable);

%% Separate EDMD/mpEDMD runs for d = 11 and d = 51

nRuns = numel(delayDimensions);
runs = repmat(struct( ...
    'dimension', [], ...
    'PsiX', [], ...
    'PsiY', [], ...
    'uniformWeight', [], ...
    'lambdaEDMD', [], ...
    'lambdaMP', [], ...
    'phiMP', []), nRuns, 1);

for runIndex = 1:nRuns
    d = delayDimensions(runIndex);
    [PsiX, PsiY] = build_delay_dictionary(observable, d);

    nRows = size(PsiX, 1);
    expectedRows = nRetainedStates - d;
    assert(nRows == expectedRows, 'Unexpected number of Hankel rows.');

    % Uniform weights are normalized over the rows of this run. A common
    % scalar weight cancels from EDMD and mpEDMD, while riggedDMD receives
    % the normalized scalar explicitly below.
    uniformWeight = 1 / nRows;

    fprintf('d = %d: %d Hankel rows; uniform row weight = %.16g\n', ...
        d, nRows, uniformWeight);

    [lambdaEDMD, lambdaMP, phiMP] = compute_edmd_mpedmd(PsiX, PsiY);

    runs(runIndex).dimension = d;
    runs(runIndex).PsiX = PsiX;
    runs(runIndex).PsiY = PsiY;
    runs(runIndex).uniformWeight = uniformWeight;
    runs(runIndex).lambdaEDMD = lambdaEDMD;
    runs(runIndex).lambdaMP = lambdaMP;
    runs(runIndex).phiMP = phiMP;

    eigenvalueFile = fullfile(outputDirectory, ...
        sprintf('logistic_eigenvalues_%d.pdf', runIndex));
    plot_eigenvalues(lambdaEDMD, lambdaMP, eigenvalueFile);
end

%% Select and plot the finite mpEDMD vector nearest lambda = -1

riggedRunIndex = find(delayDimensions == riggedDimension, 1);
assert(~isempty(riggedRunIndex), ...
    'riggedDimension must be one of the delay dimensions computed above.');

riggedRun = runs(riggedRunIndex);
[distanceToMinusOne, selectedIndex] = min(abs(riggedRun.lambdaMP + 1));
selectedEigenvalue = riggedRun.lambdaMP(selectedIndex);

fprintf(['Selected d = %d mpEDMD eigenvalue nearest -1: ' ...
    '%.16f %+.16fi; |lambda + 1| = %.6e\n'], ...
    riggedDimension, real(selectedEigenvalue), imag(selectedEigenvalue), ...
    distanceToMinusOne);

xPlot = x(1:size(riggedRun.PsiX, 1));
[xPlotSorted, sortIndex] = sort(xPlot);
selectedPhi = riggedRun.phiMP(:, selectedIndex);

plot_log_modulus(xPlotSorted, selectedPhi(sortIndex), ...
    [0.26 0.93 -14.5 -3.5], ...
    fullfile(outputDirectory, 'logistic_eigenfunction.pdf'));

%% Empirical two-cell transition frequencies on the d = 51 sample

% The partition begins at the smallest sampled attractor value rather than
% at zero. Frequencies are computed on the same 19950 one-step pairs used
% by the d = 51 dictionary and wave-packet panel.
xStar = (r - 1) / r;
nTransitionPairs = size(riggedRun.PsiX, 1);
xCurrent = x(1:nTransitionPairs);
xNext = x(2:nTransitionPairs + 1);

[transitionCounts, transitionProbabilities, flipFrequency] = ...
    two_cell_transitions(xCurrent, xNext, xStar);

fprintf('\nEmpirical two-cell partition on %d pairs:\n', nTransitionPairs);
fprintf('  L = [%.16g, %.16g), R = (%.16g, %.16g]\n', ...
    min(xCurrent), xStar, xStar, max(xCurrent));
fprintf('  Counts (rows: source L,R; columns: destination L,R):\n');
fprintf('          L          R\n');
fprintf('  L %9d  %9d\n', transitionCounts(1, 1), transitionCounts(1, 2));
fprintf('  R %9d  %9d\n', transitionCounts(2, 1), transitionCounts(2, 2));
fprintf('  Conditional frequencies P(destination | source):\n');
fprintf('          L          R\n');
fprintf('  L %9.6f  %9.6f\n', ...
    transitionProbabilities(1, 1), transitionProbabilities(1, 2));
fprintf('  R %9.6f  %9.6f\n', ...
    transitionProbabilities(2, 1), transitionProbabilities(2, 2));
fprintf('  Overall empirical flip frequency: %.6f\n\n', flipFrequency);

%% Separate riggedDMD runs for epsilon = 0.01 and epsilon = 0.75

add_riggeddmd_path();

gCoefficients = zeros(riggedDimension, 1);
gCoefficients(1) = 1;
wavePacket = [];

for epsilonIndex = 1:numel(riggedEpsilons)
    epsilonValue = riggedEpsilons(epsilonIndex);

    fprintf(['Running riggedDMD with d = %d, epsilon = %.2f, ' ...
        'theta = pi, and row weight %.16g.\n'], ...
        riggedDimension, epsilonValue, riggedRun.uniformWeight);

    [gModes, spectralDensity] = riggedDMD( ...
        riggedRun.PsiX, riggedRun.PsiY, riggedRun.uniformWeight, ...
        epsilonValue, thetaTarget, [], ...
        'order', kernelOrder, ...
        'g_coeffs', gCoefficients, ...
        'TH2', thetaGrid);

    gModes = squeeze(gModes);
    spectralDensity = real(squeeze(spectralDensity));

    densityFile = fullfile(outputDirectory, ...
        sprintf('logistic_measure_%d.pdf', epsilonIndex));

    if epsilonIndex == 1
        densityLimits = [0 pi 0 8];
    else
        densityLimits = [0 pi 0 0.6];
        % Only the (epsilon, theta) = (0.01, pi) wave packet appears in the
        % right panel of Figure 2.
        wavePacket = riggedRun.PsiX * gModes;
    end

    plot_spectral_density(thetaGrid, spectralDensity, ...
        densityLimits, densityFile);
end

assert(~isempty(wavePacket), ...
    'The epsilon = 0.01 riggedDMD wave packet was not computed.');

plot_log_modulus(xPlotSorted, wavePacket(sortIndex), ...
    [0.26 0.93 -11 0.5], ...
    fullfile(outputDirectory, 'logistic_rigged_eigenfunction.pdf'));

fprintf('Finished writing all six logistic-map panel PDFs to %s.\n', ...
    outputDirectory);

%% Local functions

function y = logistic_map(x, r)
    %LOGISTIC_MAP Evaluate the logistic map F(x) = r*x*(1-x).
    y = r * x .* (1 - x);
end

function [PsiX, PsiY] = build_delay_dictionary(observable, dimension)
    %BUILD_DELAY_DICTIONARY Construct only the required Hankel blocks.
    %
    % For m retained states and dimension d, the paired delay matrices have
    % m-d rows. Direct indexing avoids constructing an m-by-m Hankel matrix.

    nRows = numel(observable) - dimension;
    assert(nRows > 0, 'Delay dimension must be smaller than the data length.');

    rowOffsets = (1:nRows).';
    delayOffsets = 0:(dimension - 1);
    indices = rowOffsets + delayOffsets;

    PsiX = observable(indices);
    PsiY = observable(indices + 1);
end

function [lambdaEDMD, lambdaMP, phiMP] = ...
        compute_edmd_mpedmd(PsiX, PsiY)
    %COMPUTE_EDMD_MPEDMD Compute standard EDMD and QR-based mpEDMD.

    KEDMD = PsiX \ PsiY;
    lambdaEDMD = eig(KEDMD);

    [Q, R] = qr(PsiX, 'econ');
    T = (R') \ (PsiY' * Q);
    [U, ~, V] = svd(T);

    polarFactor = V * U';
    [schurVectors, schurForm] = schur(polarFactor, 'complex');

    lambdaMP = diag(schurForm);
    phiMP = PsiX * (R \ schurVectors);
end

function plot_eigenvalues(lambdaEDMD, lambdaMP, fileName)
    %PLOT_EIGENVALUES Reproduce one EDMD/mpEDMD eigenvalue panel.

    fig = figure('Color', 'w');
    hold on

    thetaCircle = linspace(0, 2*pi, 201);
    xUnit = cos(thetaCircle);
    yUnit = sin(thetaCircle);

    plot(xUnit, yUnit, 'k--', 'LineWidth', 2)
    plot(1.1*xUnit, zeros(size(yUnit)), 'k', 'LineWidth', 1)
    plot(zeros(size(xUnit)), 1.1*yUnit, 'k', 'LineWidth', 1)

    plot(real(lambdaEDMD), imag(lambdaEDMD), 'square', ...
        'LineStyle', 'none', ...
        'MarkerSize', 16, ...
        'Color', [1 69/255 79/255], ...
        'MarkerFaceColor', [1 69/255 79/255])

    plot(real(lambdaMP), imag(lambdaMP), 'diamond', ...
        'LineStyle', 'none', ...
        'MarkerSize', 16, ...
        'Color', [0 120/255 0], ...
        'MarkerFaceColor', [0 120/255 0])

    xlabel('$\mathrm{Re}(\lambda)$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    ylabel('$\mathrm{Im}(\lambda)$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    set(gca, 'FontSize', 16)
    axis([-1.1 1.1 -1.1 1.1])
    axis square
    grid on
    box on

    save_figure_pdf(fig, fileName)
end

function plot_spectral_density(thetaGrid, density, limits, fileName)
    %PLOT_SPECTRAL_DENSITY Reproduce one riggedDMD density panel.

    fig = figure('Color', 'w');
    plot(thetaGrid, density, ...
        'Color', [1 69/255 79/255], 'LineWidth', 4)

    xlabel('$\theta$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    xticks([0 pi/4 pi/2 3*pi/4 pi])
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
    set(gca, 'FontSize', 16)
    axis(limits)
    grid on
    box on

    save_figure_pdf(fig, fileName)
end

function plot_log_modulus(xSorted, valuesSorted, limits, fileName)
    %PLOT_LOG_MODULUS Reproduce an mpEDMD/riggedDMD modulus panel.

    fig = figure('Color', 'w');
    plot(xSorted, log(abs(valuesSorted)), 'k-', 'LineWidth', 1.5)

    xlabel('$x$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    ylabel('$\log|\varphi(x)|$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    set(gca, 'FontSize', 16)
    axis(limits)
    grid on
    box on

    save_figure_pdf(fig, fileName)
end

function [counts, probabilities, flipFrequency] = ...
        two_cell_transitions(xCurrent, xNext, xStar)
    %TWO_CELL_TRANSITIONS Count empirical L/R transitions.
    %
    % State 1 is L (x < xStar) and state 2 is R (x > xStar). Samples equal
    % to the boundary to roundoff are omitted rather than assigned arbitrarily.

    boundaryTolerance = 100 * eps(max(1, abs(xStar)));
    valid = abs(xCurrent - xStar) > boundaryTolerance & ...
        abs(xNext - xStar) > boundaryTolerance;

    sourceCell = 1 + (xCurrent(valid) > xStar);
    destinationCell = 1 + (xNext(valid) > xStar);

    counts = accumarray([sourceCell, destinationCell], 1, [2 2]);
    rowTotals = sum(counts, 2);
    probabilities = counts ./ rowTotals;
    flipFrequency = (counts(1, 2) + counts(2, 1)) / sum(counts(:));
end

function add_riggeddmd_path()
    %ADD_RIGGEDDMD_PATH Locate the public riggedDMD routines.

    scriptDirectory = fileparts(mfilename('fullpath'));
    candidates = { ...
        fullfile(pwd, 'main_routines'), ...
        fullfile(scriptDirectory, 'main_routines'), ...
        fullfile(scriptDirectory, '..', 'main_routines')};

    for candidateIndex = 1:numel(candidates)
        if isfolder(candidates{candidateIndex})
            addpath(candidates{candidateIndex})
            assert(exist('riggedDMD', 'file') == 2, ...
                'main_routines was found, but riggedDMD.m is missing.');
            return
        end
    end

    error(['Could not find the riggedDMD main_routines directory. ' ...
        'Place it beside this script or run the script from its parent repo.']);
end

function save_figure_pdf(figHandle, fileName)
    %SAVE_FIGURE_PDF Save a tightly sized vector PDF.

    [fileDirectory, ~, extension] = fileparts(fileName);
    if isempty(extension)
        fileName = [fileName '.pdf'];
    elseif ~strcmpi(extension, '.pdf')
        error('Figure output must have a .pdf extension.');
    end

    if ~isempty(fileDirectory) && ~isfolder(fileDirectory)
        mkdir(fileDirectory)
    end

    set(figHandle, 'PaperUnits', 'centimeters')
    set(figHandle, 'Units', 'centimeters')

    position = get(figHandle, 'Position');
    set(figHandle, 'PaperSize', [position(3) position(4)])
    set(figHandle, 'PaperPositionMode', 'manual')
    set(figHandle, 'PaperPosition', [0 0 position(3) position(4)])

    print(figHandle, fileName, '-dpdf', '-painters')
end
