% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE ROSSLER POINCARE MAP
% -------------------------------------------------------------------------
%
% This script regenerates and archives an event-located Poincare section
% for the Rossler flow, reproduces every external panel in the Rossler
% section of the paper, and records the associated numerical diagnostics.
%
% OUTPUT DATA
%   rossler_section_cache.mat       section data, provenance, diagnostics
%   rossler_section_metadata.txt    human-readable provenance/diagnostics
%
% OUTPUT FIGURES
%   rossler_2D_psec.pdf             section in the (y,z)-plane
%   rossler_1D_psec.pdf             return map y_n -> y_{n+1}
%   rossler_density.pdf             projected invariant-density histogram
%   rossler_spec_measure.pdf        d=50, epsilon=0.25 spectral density
%   rossler_eigenfunction.pdf       d=50, epsilon=0.9, theta=2*pi/3 packet
%
% The Poincare crossings are located by the event interpolation built into
% ode45. Only x=0 crossings in the positive direction are requested, and
% each retained event is checked to satisfy xdot=-y-z>0. The integration
% interval, transient, initial condition, solver, tolerances, maximum step,
% and crossing rule are stored with the section data.
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

%% Reproducibility and output parameters

regenerateSectionData = false;
cacheFile = 'rossler_section_cache.mat';
metadataFile = 'rossler_section_metadata.txt';
outputDirectory = '.';

if ~isfolder(outputDirectory)
    mkdir(outputDirectory)
end

%% Rossler and integration parameters

config = struct();
config.a = 0.1;
config.b = 0.1;
config.c = 18;
config.initialCondition = [0.1; -8; 0.3];
config.integrationInterval = [0, 100000];
config.transientDiscardTime = 100;
config.solver = 'ode45';
config.relativeTolerance = 1e-12;
config.absoluteTolerance = [1e-12, 1e-12, 1e-12];
config.maximumStep = 0.1;
config.sectionCoordinate = 'x';
config.sectionValue = 0;
config.eventDirection = +1;
config.crossingRule = [ ...
    'Event-located x=0 crossing with positive direction; retain only ' ...
    'events satisfying xdot=-y-z>0 after the discarded transient.'];

%% riggedDMD parameters

delayDimension = 50;
spectralEpsilon = 0.25;
wavePacketEpsilon = 0.9;
wavePacketTheta = 2*pi/3;
thetaGrid = -pi:0.01:pi;
kernelOrder = 2;

%% Generate or load one authoritative section cache

if regenerateSectionData || ~isfile(cacheFile)
    cache = generate_section_cache(config);
    save(cacheFile, 'cache', '-v7')
else
    loaded = load(cacheFile, 'cache');
    assert(isfield(loaded, 'cache'), ...
        'The cache file does not contain a variable named cache.');
    cache = loaded.cache;
    validate_cache(cache, config);
end

poincareXYZ = cache.retainedCrossingsXYZ;
crossingTimes = cache.retainedCrossingTimes;
yData = poincareXYZ(:, 2);
zData = poincareXYZ(:, 3);
N = numel(yData);

assert(N > delayDimension, ...
    'The retained section must contain more points than the delay dimension.');

fprintf('Retained Poincare intersections N = %d\n', N);
fprintf('Discarded events before t = %.16g: %d\n', ...
    config.transientDiscardTime, cache.discardedEventCount);
fprintf('Retained crossing-time interval: [%.16g, %.16g]\n', ...
    crossingTimes(1), crossingTimes(end));
fprintf('Maximum retained |x| at an event: %.6e\n', ...
    max(abs(poincareXYZ(:, 1))));
fprintf('Minimum retained xdot: %.16g\n', ...
    min(cache.retainedCrossingSpeeds));

%% Re-estimate the fixed point and its nontrivial preimage

yCurrent = yData(1:end-1);
yNext = yData(2:end);

[fixedPointY, fixedPointInfo] = estimate_fixed_point(yCurrent, yNext);
[fixedPointPreimageY, preimageInfo] = ...
    estimate_right_preimage(yCurrent, yNext, fixedPointY);

fprintf('Estimated fixed point y_* = %.10f\n', fixedPointY);
fprintf('Estimated nontrivial preimage y_** = %.10f\n', ...
    fixedPointPreimageY);

%% Re-evaluate the empirical three-cell transition support

[transitionCounts, transitionProbabilities, transitionSupport] = ...
    three_cell_transitions(yData, fixedPointY, fixedPointPreimageY);

fprintf('\nThree-cell order: (I1,I2,I3)\n');
fprintf('  I1: y > y_**, I2: y_* < y <= y_**, I3: y <= y_*\n');
fprintf('Empirical transition counts:\n');
disp(transitionCounts)
fprintf('Empirical conditional transition probabilities:\n');
disp(transitionProbabilities)
fprintf('Empirical binary transition support (count > 0):\n');
disp(transitionSupport)

%% Re-evaluate the projected histogram and nearby peaks

[histogramDensity, histogramEdges] = histcounts( ...
    yData, 50, 'Normalization', 'pdf');
histogramCenters = (histogramEdges(1:end-1) + histogramEdges(2:end)) / 2;

[histogramPeakNearFixedPoint, histogramValueNearFixedPoint] = ...
    nearby_maximum(histogramCenters, histogramDensity, fixedPointY, 1.0);
[histogramPeakNearPreimage, histogramValueNearPreimage] = ...
    nearby_maximum(histogramCenters, histogramDensity, ...
    fixedPointPreimageY, 1.0);

%% Plot the Poincare section in the (y,z)-plane

fig = figure('Color', 'w');
plot(yData, zData, 'k.', 'MarkerSize', 10)
set(gca, 'FontSize', 16)
xlabel('$y$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
ylabel('$z$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-28 -11 0.00505 0.0055])

save_figure_pdf(fig, ...
    fullfile(outputDirectory, 'rossler_2D_psec.pdf'))

%% Plot the one-dimensional return map in y

fig = figure('Color', 'w');
hold on

plot(yCurrent, yCurrent, '.', ...
    'Color', [0.5 0.5 0.5], 'MarkerSize', 5)
plot(yCurrent, yNext, 'k.', 'MarkerSize', 10)

plot(fixedPointY, fixedPointY, 'square', ...
    'LineStyle', 'none', ...
    'MarkerSize', 20, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

plot(fixedPointPreimageY, fixedPointY, 'diamond', ...
    'LineStyle', 'none', ...
    'MarkerSize', 20, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

set(gca, 'FontSize', 16)
xlabel('$y_n$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
ylabel('$y_{n+1}$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-27 -11.5 -27 -11.5])

save_figure_pdf(fig, ...
    fullfile(outputDirectory, 'rossler_1D_psec.pdf'))

%% Plot the projected invariant-density histogram

fig = figure('Color', 'w');
histogram(yData, histogramEdges, ...
    'Normalization', 'pdf', ...
    'FaceColor', [36/255 122/255 254/255])

set(gca, 'FontSize', 16)
xlabel('$y$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
ylabel('Density', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
axis([-27 -11.5 0 0.25])
box on

save_figure_pdf(fig, ...
    fullfile(outputDirectory, 'rossler_density.pdf'))

%% Build d=50 delay-coordinate matrices from the same cache

centeredY = yData - mean(yData);
normalizedY = centeredY / std(centeredY);

[PsiX, PsiY] = build_delay_dictionary(normalizedY, delayDimension);
nSnapshotPairs = size(PsiX, 1);

assert(nSnapshotPairs == N - delayDimension, ...
    'The d=50 delay matrices must have N-d rows.');

% Use an explicit uniform quadrature vector so that the weights sum to one
% over exactly the N-d paired Hankel snapshots.
uniformWeights = ones(nSnapshotPairs, 1) / nSnapshotPairs;
assert(abs(sum(uniformWeights) - 1) < 1e-12, ...
    'Uniform quadrature weights are not normalized.');

fprintf('\nd = %d delay-coordinate dictionary\n', delayDimension);
fprintf('Hankel snapshot pairs N-d = %d\n', nSnapshotPairs);
fprintf('Uniform quadrature weight = %.16g; sum(weights) = %.16g\n', ...
    uniformWeights(1), sum(uniformWeights));

%% Locate riggedDMD and define the observable coefficients

add_riggeddmd_path();

gCoefficients = zeros(delayDimension, 1);
gCoefficients(1) = 1;

%% Separate riggedDMD run for the spectral density: epsilon=0.25

fprintf('Running spectral-density riggedDMD: d=%d, epsilon=%.2f.\n', ...
    delayDimension, spectralEpsilon);

[~, spectralDensity] = riggedDMD( ...
    PsiX, PsiY, uniformWeights, spectralEpsilon, wavePacketTheta, [], ...
    'order', kernelOrder, ...
    'g_coeffs', gCoefficients, ...
    'TH2', thetaGrid);

spectralDensity = real(squeeze(spectralDensity));

positiveThetaMask = thetaGrid >= 0 & thetaGrid <= pi;
positiveTheta = thetaGrid(positiveThetaMask);
positiveDensity = spectralDensity(positiveThetaMask);
[spectralPeakValue, spectralPeakIndex] = max(positiveDensity);
spectralPeakTheta = positiveTheta(spectralPeakIndex);

fprintf('Spectral-density peak on [0,pi]: theta = %.10f = %.10f*pi\n', ...
    spectralPeakTheta, spectralPeakTheta/pi);
fprintf('Spectral-density peak value: %.10f\n', spectralPeakValue);

fig = figure('Color', 'w');
plot(thetaGrid, spectralDensity, ...
    'Color', [1 69/255 79/255], 'LineWidth', 4)

xlabel('$\theta$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
set(gca, 'FontSize', 16)
axis([0 pi 0 0.35])
grid on
box on

save_figure_pdf(fig, ...
    fullfile(outputDirectory, 'rossler_spec_measure.pdf'))

%% Separate riggedDMD run for the wave packet: epsilon=0.9, theta=2*pi/3

fprintf(['Running wave-packet riggedDMD: d=%d, epsilon=%.2f, ' ...
    'theta=%.10f.\n'], ...
    delayDimension, wavePacketEpsilon, wavePacketTheta);

[wavePacketCoefficients, ~] = riggedDMD( ...
    PsiX, PsiY, uniformWeights, wavePacketEpsilon, ...
    wavePacketTheta, [], ...
    'order', kernelOrder, ...
    'g_coeffs', gCoefficients, ...
    'TH2', []);

wavePacketCoefficients = squeeze(wavePacketCoefficients);
wavePacket = PsiX * wavePacketCoefficients;
wavePacketLogModulus = log(abs(wavePacket));

% Plot against the original, unnormalized y-coordinate corresponding to
% the first entry in each delay row.
yPlot = yData(1:nSnapshotPairs);
[yPlotSorted, sortIndex] = sort(yPlot);
wavePacketLogModulusSorted = wavePacketLogModulus(sortIndex);

[troughNearFixedPoint, troughValueNearFixedPoint] = nearby_minimum( ...
    yPlot, wavePacketLogModulus, fixedPointY, 1.0);
[troughNearPreimage, troughValueNearPreimage] = nearby_minimum( ...
    yPlot, wavePacketLogModulus, fixedPointPreimageY, 1.0);

fprintf('Wave-packet trough near y_*: y = %.10f, log-modulus = %.10f\n', ...
    troughNearFixedPoint, troughValueNearFixedPoint);
fprintf('Wave-packet trough near y_**: y = %.10f, log-modulus = %.10f\n', ...
    troughNearPreimage, troughValueNearPreimage);

fig = figure('Color', 'w');
hold on

xline(fixedPointY, '--', ...
    'Color', [1 69/255 79/255], 'LineWidth', 2)
xline(fixedPointPreimageY, '--', ...
    'Color', [0 120/255 0], 'LineWidth', 2)

plot(yPlotSorted, wavePacketLogModulusSorted, ...
    'k-', 'LineWidth', 1.5)

% This abscissa is the unnormalized section coordinate y, not x.
xlabel('$y$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\log|\varphi(y)|$', 'Interpreter', 'latex', ...
    'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
axis([-27 -11.5 -7 0])
box on
grid on

save_figure_pdf(fig, ...
    fullfile(outputDirectory, 'rossler_eigenfunction.pdf'))

%% Archive all derived values and write a readable report

analysis = struct();
analysis.generatedOn = datestr(now, 30);
analysis.N = N;
analysis.delayDimension = delayDimension;
analysis.nSnapshotPairs = nSnapshotPairs;
analysis.uniformWeight = uniformWeights(1);
analysis.uniformWeightSum = sum(uniformWeights);
analysis.fixedPointY = fixedPointY;
analysis.fixedPointFit = fixedPointInfo;
analysis.fixedPointPreimageY = fixedPointPreimageY;
analysis.preimageFit = preimageInfo;
analysis.transitionCellOrder = {'I1', 'I2', 'I3'};
analysis.transitionCounts = transitionCounts;
analysis.transitionProbabilities = transitionProbabilities;
analysis.transitionSupport = transitionSupport;
analysis.histogramBinCount = 50;
analysis.histogramPeakNearFixedPoint = histogramPeakNearFixedPoint;
analysis.histogramValueNearFixedPoint = histogramValueNearFixedPoint;
analysis.histogramPeakNearPreimage = histogramPeakNearPreimage;
analysis.histogramValueNearPreimage = histogramValueNearPreimage;
analysis.spectralEpsilon = spectralEpsilon;
analysis.spectralPeakTheta = spectralPeakTheta;
analysis.spectralPeakThetaOverPi = spectralPeakTheta/pi;
analysis.spectralPeakValue = spectralPeakValue;
analysis.wavePacketEpsilon = wavePacketEpsilon;
analysis.wavePacketTheta = wavePacketTheta;
analysis.wavePacketTroughNearFixedPoint = troughNearFixedPoint;
analysis.wavePacketTroughValueNearFixedPoint = troughValueNearFixedPoint;
analysis.wavePacketTroughNearPreimage = troughNearPreimage;
analysis.wavePacketTroughValueNearPreimage = troughValueNearPreimage;
analysis.riggedDMDImplementation = which('riggedDMD');
analysis.mpEDMDImplementation = which('mpEDMDqr');

cache.analysis = analysis;
save(cacheFile, 'cache', '-v7')
write_metadata_report(metadataFile, cache);

fprintf('\nFinished writing the Rossler cache, diagnostics, and five PDFs.\n');

%% Local functions

function cache = generate_section_cache(config)
    %GENERATE_SECTION_CACHE Integrate once and retain event-located crossings.

    eventFunction = @(t, state) positive_x_section_event(t, state);
    odeOptions = odeset( ...
        'RelTol', config.relativeTolerance, ...
        'AbsTol', config.absoluteTolerance, ...
        'MaxStep', config.maximumStep, ...
        'Events', eventFunction);

    fprintf(['Generating event-located section data on [%.16g, %.16g] ' ...
        'with %s.\n'], ...
        config.integrationInterval(1), config.integrationInterval(2), ...
        config.solver);

    [~, ~, eventTimes, eventStates, eventIndices] = ode45( ...
        @(t, state) rossler_rhs(state, config.a, config.b, config.c), ...
        config.integrationInterval, config.initialCondition, odeOptions);

    eventSpeeds = -eventStates(:, 2) - eventStates(:, 3);
    positiveCrossing = eventSpeeds > 0;
    afterTransient = eventTimes >= config.transientDiscardTime;
    retainedMask = positiveCrossing & afterTransient;

    assert(all(eventIndices == 1), ...
        'Unexpected event index returned by ode45.');
    assert(any(retainedMask), 'No Poincare crossings survived the filters.');

    cache = struct();
    cache.schemaVersion = 1;
    cache.createdOn = datestr(now, 30);
    cache.system = 'Rossler';
    cache.parameters = struct('a', config.a, 'b', config.b, 'c', config.c);
    cache.initialCondition = config.initialCondition;
    cache.integrationInterval = config.integrationInterval;
    cache.transientDiscardTime = config.transientDiscardTime;
    cache.solver = config.solver;
    cache.relativeTolerance = config.relativeTolerance;
    cache.absoluteTolerance = config.absoluteTolerance;
    cache.maximumStep = config.maximumStep;
    cache.sectionCoordinate = config.sectionCoordinate;
    cache.sectionValue = config.sectionValue;
    cache.eventDirection = config.eventDirection;
    cache.crossingRule = config.crossingRule;
    cache.eventLocationRule = ...
        'Root located by the continuous event interpolation used by ode45.';
    cache.rawEventTimes = eventTimes;
    cache.rawEventStatesXYZ = eventStates;
    cache.rawEventSpeeds = eventSpeeds;
    cache.retainedEventMask = retainedMask;
    cache.discardedEventCount = sum(~retainedMask);
    cache.retainedCrossingTimes = eventTimes(retainedMask);
    cache.retainedCrossingsXYZ = eventStates(retainedMask, :);
    cache.retainedCrossingSpeeds = eventSpeeds(retainedMask);
    cache.N = sum(retainedMask);
end

function validate_cache(cache, config)
    %VALIDATE_CACHE Reject an undocumented or parameter-mismatched cache.

    requiredFields = { ...
        'schemaVersion', 'parameters', 'initialCondition', ...
        'integrationInterval', 'transientDiscardTime', 'solver', ...
        'relativeTolerance', 'absoluteTolerance', 'maximumStep', ...
        'crossingRule', 'retainedCrossingTimes', ...
        'retainedCrossingsXYZ', 'retainedCrossingSpeeds', ...
        'discardedEventCount', 'N'};

    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(cache, requiredFields{fieldIndex}), ...
            ['Cache lacks required provenance field: ' ...
            requiredFields{fieldIndex}]);
    end

    matches = ...
        cache.parameters.a == config.a && ...
        cache.parameters.b == config.b && ...
        cache.parameters.c == config.c && ...
        isequal(cache.initialCondition(:), config.initialCondition(:)) && ...
        isequal(cache.integrationInterval, config.integrationInterval) && ...
        cache.transientDiscardTime == config.transientDiscardTime && ...
        strcmp(cache.solver, config.solver) && ...
        cache.relativeTolerance == config.relativeTolerance && ...
        isequal(cache.absoluteTolerance, config.absoluteTolerance) && ...
        cache.maximumStep == config.maximumStep;

    assert(matches, [ ...
        'The existing section cache does not match the requested settings. ' ...
        'Set regenerateSectionData=true to rebuild it.']);

    assert(cache.N == size(cache.retainedCrossingsXYZ, 1), ...
        'The cached N does not match the retained section data.');
    assert(all(cache.retainedCrossingSpeeds > 0), ...
        'The cache contains a crossing that violates xdot>0.');
end

function [value, isTerminal, direction] = ...
        positive_x_section_event(~, state)
    %POSITIVE_X_SECTION_EVENT Locate x=0 crossings in the +x direction.
    value = state(1);
    isTerminal = 0;
    direction = +1;
end

function derivative = rossler_rhs(state, a, b, c)
    %ROSSLER_RHS Right-hand side of the Rossler system.
    derivative = [ ...
        -state(2) - state(3); ...
         state(1) + a*state(2); ...
         b + state(3)*(state(1) - c)];
end

function [fixedPoint, info] = estimate_fixed_point(yCurrent, yNext)
    %ESTIMATE_FIXED_POINT Estimate F(y)=y by iterated local linear fits.
    residual = yNext - yCurrent;
    eligible = true(size(yCurrent));
    [fixedPoint, info] = refine_zero( ...
        yCurrent, residual, eligible, [0.5, 0.25, 0.1]);
end

function [preimage, info] = ...
        estimate_right_preimage(yCurrent, yNext, fixedPoint)
    %ESTIMATE_RIGHT_PREIMAGE Estimate the nontrivial right preimage of y_*.
    residual = yNext - fixedPoint;
    eligible = yCurrent > fixedPoint + 1;
    [preimage, info] = refine_zero( ...
        yCurrent, residual, eligible, [0.75, 0.35, 0.15]);
end

function [rootEstimate, info] = ...
        refine_zero(abscissa, residual, eligible, halfWidths)
    %REFINE_ZERO Refine a sampled zero using centered local linear fits.

    candidateResidual = abs(residual);
    candidateResidual(~eligible) = inf;
    [minimumResidual, candidateIndex] = min(candidateResidual);
    assert(isfinite(minimumResidual), 'No eligible root candidate was found.');

    rootEstimate = abscissa(candidateIndex);
    fitCounts = zeros(size(halfWidths));
    slopes = zeros(size(halfWidths));

    for widthIndex = 1:numel(halfWidths)
        width = halfWidths(widthIndex);
        localMask = eligible & abs(abscissa - rootEstimate) <= width;
        fitCounts(widthIndex) = sum(localMask);
        assert(fitCounts(widthIndex) >= 10, ...
            'Too few samples for a local return-map fit.');

        centeredAbscissa = abscissa(localMask) - rootEstimate;
        fitCoefficients = polyfit(centeredAbscissa, residual(localMask), 1);
        slopes(widthIndex) = fitCoefficients(1);
        assert(abs(slopes(widthIndex)) > sqrt(eps), ...
            'Local return-map residual has an unresolved slope.');

        rootEstimate = rootEstimate - ...
            fitCoefficients(2) / fitCoefficients(1);
    end

    info = struct();
    info.initialCandidate = abscissa(candidateIndex);
    info.initialAbsoluteResidual = minimumResidual;
    info.halfWidths = halfWidths;
    info.fitCounts = fitCounts;
    info.localSlopes = slopes;
end

function [counts, probabilities, support] = ...
        three_cell_transitions(yData, fixedPoint, preimage)
    %THREE_CELL_TRANSITIONS Empirical transitions in order (I1,I2,I3).

    labels = zeros(size(yData));
    labels(yData > preimage) = 1;
    labels(yData > fixedPoint & yData <= preimage) = 2;
    labels(yData <= fixedPoint) = 3;

    assert(all(labels > 0), 'A section point was not assigned to a cell.');

    source = labels(1:end-1);
    destination = labels(2:end);
    counts = accumarray([source, destination], 1, [3 3]);

    rowTotals = sum(counts, 2);
    probabilities = bsxfun(@rdivide, counts, rowTotals);
    support = counts > 0;
end

function [location, value] = ...
        nearby_maximum(abscissa, ordinate, target, halfWidth)
    %NEARBY_MAXIMUM Find a sampled local maximum near a target.
    localMask = abs(abscissa - target) <= halfWidth & isfinite(ordinate);
    assert(any(localMask), 'No finite sample lies in the peak-search window.');
    localIndices = find(localMask);
    [value, relativeIndex] = max(ordinate(localMask));
    absoluteIndex = localIndices(relativeIndex);
    location = abscissa(absoluteIndex);
end

function [location, value] = ...
        nearby_minimum(abscissa, ordinate, target, halfWidth)
    %NEARBY_MINIMUM Find a sampled local minimum near a target.
    localMask = abs(abscissa - target) <= halfWidth & isfinite(ordinate);
    assert(any(localMask), 'No finite sample lies in the trough-search window.');
    localIndices = find(localMask);
    [value, relativeIndex] = min(ordinate(localMask));
    absoluteIndex = localIndices(relativeIndex);
    location = abscissa(absoluteIndex);
end

function [PsiX, PsiY] = build_delay_dictionary(observable, dimension)
    %BUILD_DELAY_DICTIONARY Construct only the required Hankel blocks.

    nRows = numel(observable) - dimension;
    assert(nRows > 0, 'Delay dimension must be smaller than the data length.');

    rowOffsets = (1:nRows).';
    delayOffsets = 0:(dimension - 1);
    indices = rowOffsets + delayOffsets;

    PsiX = observable(indices);
    PsiY = observable(indices + 1);
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
            break
        end
    end

    assert(exist('riggedDMD', 'file') == 2, ...
        'Could not find riggedDMD.m in main_routines or on the MATLAB path.');
    assert(exist('mpEDMDqr', 'file') == 2, ...
        'Could not find the mpEDMDqr.m dependency on the MATLAB path.');
end

function write_metadata_report(fileName, cache)
    %WRITE_METADATA_REPORT Write provenance and diagnostics as plain text.

    fileID = fopen(fileName, 'w');
    assert(fileID >= 0, 'Could not open the metadata report for writing.');
    closeFile = onCleanup(@() fclose(fileID)); 

    analysis = cache.analysis;

    fprintf(fileID, 'Rossler Poincare-section archive\n');
    fprintf(fileID, 'Created: %s\n', cache.createdOn);
    fprintf(fileID, 'Analyzed: %s\n\n', analysis.generatedOn);

    fprintf(fileID, 'System parameters: a=%.16g, b=%.16g, c=%.16g\n', ...
        cache.parameters.a, cache.parameters.b, cache.parameters.c);
    fprintf(fileID, 'Initial condition: [%.16g, %.16g, %.16g]\n', ...
        cache.initialCondition(1), cache.initialCondition(2), ...
        cache.initialCondition(3));
    fprintf(fileID, 'Integration interval: [%.16g, %.16g]\n', ...
        cache.integrationInterval(1), cache.integrationInterval(2));
    fprintf(fileID, 'Discarded transient: t < %.16g\n', ...
        cache.transientDiscardTime);
    fprintf(fileID, 'Discarded events: %d\n', cache.discardedEventCount);
    fprintf(fileID, 'Solver: %s\n', cache.solver);
    fprintf(fileID, 'Relative tolerance: %.16g\n', ...
        cache.relativeTolerance);
    fprintf(fileID, 'Absolute tolerances: [%.16g, %.16g, %.16g]\n', ...
        cache.absoluteTolerance(1), cache.absoluteTolerance(2), ...
        cache.absoluteTolerance(3));
    fprintf(fileID, 'Maximum solver step: %.16g\n', cache.maximumStep);
    fprintf(fileID, 'Crossing rule: %s\n', cache.crossingRule);
    fprintf(fileID, 'Event-location rule: %s\n\n', cache.eventLocationRule);

    fprintf(fileID, 'Retained section points N: %d\n', analysis.N);
    fprintf(fileID, 'Delay dimension d: %d\n', analysis.delayDimension);
    fprintf(fileID, 'Hankel snapshot pairs N-d: %d\n', ...
        analysis.nSnapshotPairs);
    fprintf(fileID, 'Uniform weight: %.16g\n', analysis.uniformWeight);
    fprintf(fileID, 'Sum of weights: %.16g\n\n', ...
        analysis.uniformWeightSum);

    fprintf(fileID, 'Estimated y_*: %.10f\n', analysis.fixedPointY);
    fprintf(fileID, 'Estimated y_**: %.10f\n', ...
        analysis.fixedPointPreimageY);
    fprintf(fileID, 'Histogram peak near y_*: %.10f (density %.10f)\n', ...
        analysis.histogramPeakNearFixedPoint, ...
        analysis.histogramValueNearFixedPoint);
    fprintf(fileID, 'Histogram peak near y_**: %.10f (density %.10f)\n\n', ...
        analysis.histogramPeakNearPreimage, ...
        analysis.histogramValueNearPreimage);

    fprintf(fileID, 'Transition-cell order: (I1,I2,I3)\n');
    fprintf(fileID, 'Transition counts:\n%s\n', ...
        matrix_as_text(analysis.transitionCounts, '%12d'));
    fprintf(fileID, 'Conditional transition probabilities:\n%s\n', ...
        matrix_as_text(analysis.transitionProbabilities, '%12.8f'));
    fprintf(fileID, 'Binary transition support:\n%s\n', ...
        matrix_as_text(double(analysis.transitionSupport), '%12d'));

    fprintf(fileID, 'Spectral epsilon: %.16g\n', analysis.spectralEpsilon);
    fprintf(fileID, 'Spectral peak theta: %.10f = %.10f*pi\n', ...
        analysis.spectralPeakTheta, analysis.spectralPeakThetaOverPi);
    fprintf(fileID, 'Spectral peak value: %.10f\n\n', ...
        analysis.spectralPeakValue);

    fprintf(fileID, 'Wave-packet epsilon: %.16g\n', ...
        analysis.wavePacketEpsilon);
    fprintf(fileID, 'Wave-packet theta: %.10f\n', ...
        analysis.wavePacketTheta);
    fprintf(fileID, 'Trough near y_*: %.10f (log-modulus %.10f)\n', ...
        analysis.wavePacketTroughNearFixedPoint, ...
        analysis.wavePacketTroughValueNearFixedPoint);
    fprintf(fileID, 'Trough near y_**: %.10f (log-modulus %.10f)\n\n', ...
        analysis.wavePacketTroughNearPreimage, ...
        analysis.wavePacketTroughValueNearPreimage);

    fprintf(fileID, 'riggedDMD implementation: %s\n', ...
        analysis.riggedDMDImplementation);
    fprintf(fileID, 'mpEDMD implementation: %s\n', ...
        analysis.mpEDMDImplementation);
end

function output = matrix_as_text(matrix, format)
    %MATRIX_AS_TEXT Convert a small numeric matrix to aligned text.

    output = '';
    for rowIndex = 1:size(matrix, 1)
        for columnIndex = 1:size(matrix, 2)
            output = [output sprintf(format, ...
                matrix(rowIndex, columnIndex))]; 
        end
        if rowIndex < size(matrix, 1)
            output = [output sprintf('\n')];
        end
    end
end

function save_figure_pdf(figHandle, fileName)
    %SAVE_FIGURE_PDF Save a figure to a tightly sized vector PDF.

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
