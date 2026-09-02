% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE FORCED DUFFING POINCARE MAP: REGIME II
% -------------------------------------------------------------------------
%
% This script generates or reuses an event-located stroboscopic section
% for the forced Duffing oscillator and analyzes its dominant seven-region
% organization with rigged Dynamic Mode Decomposition (riggedDMD).
%
% The autonomous phase formulation is
%
%   xdot     = v,
%   vdot     = -delta*v - alpha*x - beta*x^3 + gamma*cos(phi),
%   phidot   = omega.
%
% Stroboscopic intersections are roots of sin(phi)=0 with positive
% direction, equivalently phi=0 mod 2*pi. MATLAB's continuous event
% interpolation locates each crossing. The trajectory is integrated in
% blocks to bound memory use, with the terminal state of each block passed
% unchanged to the next block.
%
% The analysis uses the standardized position observable and d=20 delays.
% It performs the smoothing sweep reported in the paper, orients each
% finite wave packet to have mean phase advance near +6*pi/7, and computes
% the following diagnostics:
%
%   - the local spectral maximum near 6*pi/7;
%   - target mismatch and phase-residual concentration;
%   - reliable one-step agreement with R_j -> R_{j+3} mod 7;
%   - the unfiltered seven-step same-region rate;
%   - enrichment of transport errors in the lowest retained modulus
%     decile; and
%   - spectral robustness across x, v, and x+v observables.
%
% The endpoint-wise reliability filter retains a transition only when the
% wave-packet modulus at both endpoints exceeds its fifth percentile. The
% quadrature weight is normalized over the N-d Hankel snapshot pairs.
%
% REQUIRED CODE
%   main_routines/riggedDMD.m and its dependencies
%
% OUTPUT DATA
%   duffing_section_II_cache.mat
%       Event-located section data, numerical provenance, spectra,
%       wave packets, transport diagnostics, and robustness results
%   duffing_section_II_metadata.txt
%       Human-readable generation and analysis metadata
%
% OUTPUT FIGURES
%   duffing_regime2_spectral_density.pdf
%       Spectral density at epsilon=0.25, displaying the three peaks
%   duffing_regime2_phase_partition.pdf
%       Seven-region phase partition at epsilon=0.10
%   duffing_regime2_one_step_transport.pdf
%       The images F(R_j) at epsilon=0.10, coloured by source region
%
% All output files are written beside this script. No output directory is
% created.
%
% The riggedDMD implementation is available from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% Paths and reproducibility controls

scriptDirectory = fileparts(mfilename('fullpath'));
if isempty(scriptDirectory)
    scriptDirectory = pwd;
end

outputDirectory = scriptDirectory;
cacheFile = fullfile(outputDirectory, 'duffing_section_II_cache.mat');
metadataFile = fullfile( ...
    outputDirectory, 'duffing_section_II_metadata.txt');
routinesDirectory = fullfile(scriptDirectory, 'main_routines');

forceRegenerateSectionData = false;
saveFigures = true;
runObservableRobustness = true;

%% Duffing model and section-generation parameters

config = struct();
config.schemaVersion = 2;
config.system = 'Forced Duffing oscillator, Regime II';
config.alpha = -1;
config.beta = 0.25;
config.delta = 0.1;
config.omega = 2;
config.gamma = 2.5;
config.initialCondition = [1; 0; 0];
config.outputStep = 1e-3;
config.numberOfOriginalOutputTimes = 1e8;
config.integrationInterval = [0, ...
    (config.numberOfOriginalOutputTimes - 1)*config.outputStep];
config.solver = 'ode45';
config.relativeTolerance = 1e-12;
config.absoluteTolerance = 1e-12*ones(1, 3);
config.maximumStep = config.outputStep;
config.refine = 1;
config.integrationBlockLength = 500;
config.discardedCrossings = 0;
config.forcingPhase = 0;
config.eventDirection = +1;
config.phasePeriod = 2*pi;
config.forcingPeriod = 2*pi/config.omega;
config.crossingRule = [ ...
    'Event-located root of sin(phi)=0 with positive direction, ' ...
    'excluding the initial point; equivalently phi=0 mod 2*pi.'];
config.eventLocationRule = [ ...
    'Root located by the continuous event interpolation used by ode45.'];

%% riggedDMD and transport parameters

delayDimension = 20;
epsilonSweep = [0.05, 0.10, 0.15, 0.20, 0.25];
spectralFigureEpsilon = 0.25;
phaseFigureEpsilon = 0.10;
kernelOrder = 2;
thetaTarget = 6*pi/7;
thetaGrid = linspace(-pi, pi, 1601);
numberOfRegions = 7;

reliabilityQuantile = 0.05;
lowModulusQuantile = 0.10;

%% Generate or load the documented stroboscopic section

[x, generationMetadata, generatedNow] = load_or_generate_section_data( ...
    cacheFile, config, forceRegenerateSectionData);

N = size(x, 1);
initialRelativePhase = ...
    config.initialCondition(3) - config.forcingPhase;
finalRelativePhase = initialRelativePhase + ...
    config.omega*diff(config.integrationInterval);
expectedRawCrossings = floor( ...
    finalRelativePhase/config.phasePeriod) - ...
    floor(initialRelativePhase/config.phasePeriod);
expectedRetainedPoints = ...
    expectedRawCrossings - config.discardedCrossings;

assert(N == expectedRetainedPoints, ...
    'Expected %d retained samples, but the cache contains %d.', ...
    expectedRetainedPoints, N)
assert(N == 31830, ...
    'The documented Regime II computation must retain N=31830 samples.')

fprintf('\nDuffing Regime II section data\n')
fprintf('------------------------------\n')
fprintf('Retained stroboscopic intersections N = %d\n', N)
fprintf('Discarded initial intersections = %d\n', ...
    generationMetadata.discardedCrossings)
fprintf('Integration interval = [%.16g, %.16g]\n', ...
    generationMetadata.integrationInterval)
fprintf('Retained time interval = [%.16g, %.16g]\n', ...
    generationMetadata.retainedEventTimes(1), ...
    generationMetadata.retainedEventTimes(end))
fprintf('Maximum event phase residual = %.6e rad\n', ...
    generationMetadata.maximumEventPhaseResidual)
if generatedNow
    fprintf('Saved the regenerated event-located cache to %s\n', cacheFile)
else
    fprintf('Loaded the documented event-located cache from %s\n', cacheFile)
end

%% Locate the riggedDMD routines

if isfolder(routinesDirectory)
    addpath(routinesDirectory)
end
assert(exist('riggedDMD', 'file') == 2, ...
    ['Could not find riggedDMD.m. Add the riggedDMD main_routines ' ...
     'directory to the MATLAB path.'])

%% Build the position delay-coordinate dictionary

[positionObservable, positionMean, positionScale] = ...
    standardize_observable(x(:, 1));
[PsiX, PsiY] = build_delay_dictionary( ...
    positionObservable, delayDimension);

nSnapshotPairs = size(PsiX, 1);
sectionData = x(1:nSnapshotPairs, :);

assert(nSnapshotPairs == N - delayDimension, ...
    'The delay matrices must contain N-d paired snapshots.')
assert(nSnapshotPairs == 31810, ...
    'The documented Regime II computation must use M=31810 snapshots.')

uniformWeights = ones(nSnapshotPairs, 1)/nSnapshotPairs;
uniformWeightSum = sum(uniformWeights);
assert(abs(uniformWeightSum - 1) < 1e-12, ...
    'The uniform quadrature weights do not sum to one.')

gCoefficients = zeros(delayDimension, 1);
gCoefficients(1) = 1;

fprintf('\nd = %d delay-coordinate dictionary\n', delayDimension)
fprintf('Hankel snapshot pairs N-d = %d\n', nSnapshotPairs)
fprintf('Uniform quadrature weight = %.16g; sum(weights) = %.16g\n', ...
    uniformWeights(1), uniformWeightSum)

%% Sweep the smoothing parameter

numberOfEpsilonValues = numel(epsilonSweep);
spectralDensitySweep = nan(numel(thetaGrid), numberOfEpsilonValues);
peakLocation = nan(numberOfEpsilonValues, 1);
peakHeight = nan(numberOfEpsilonValues, 1);
meanPhaseAdvance = nan(numberOfEpsilonValues, 1);
targetMismatch = nan(numberOfEpsilonValues, 1);
phaseConcentration = nan(numberOfEpsilonValues, 1);
medianAbsoluteResidual = nan(numberOfEpsilonValues, 1);
transitionAccuracy = nan(numberOfEpsilonValues, 1);
sevenStepReturn = nan(numberOfEpsilonValues, 1);
baselineTransportError = nan(numberOfEpsilonValues, 1);
lowModulusTransportError = nan(numberOfEpsilonValues, 1);
transportErrorEnrichment = nan(numberOfEpsilonValues, 1);
metricData = cell(numberOfEpsilonValues, 1);

for epsilonIndex = 1:numberOfEpsilonValues
    epsilon = epsilonSweep(epsilonIndex);
    fprintf('Running riggedDMD with epsilon = %.2f.\n', epsilon)

    [gMode, spectralDensity] = run_rigged_dmd( ...
        PsiX, PsiY, uniformWeights, epsilon, thetaTarget, ...
        kernelOrder, gCoefficients, thetaGrid);

    rawWavePacket = PsiX*gMode;
    metrics = transport_metrics( ...
        rawWavePacket, thetaTarget, reliabilityQuantile, ...
        lowModulusQuantile, numberOfRegions);

    spectralDensitySweep(:, epsilonIndex) = real(spectralDensity(:));
    [peakLocation(epsilonIndex), peakHeight(epsilonIndex)] = ...
        local_maximum( ...
            thetaGrid, spectralDensitySweep(:, epsilonIndex), ...
            thetaTarget, 0.06*pi);
    meanPhaseAdvance(epsilonIndex) = metrics.meanPhaseAdvance;
    targetMismatch(epsilonIndex) = metrics.targetMismatch;
    phaseConcentration(epsilonIndex) = metrics.phaseConcentration;
    medianAbsoluteResidual(epsilonIndex) = ...
        metrics.medianAbsoluteResidual;
    transitionAccuracy(epsilonIndex) = metrics.transitionAccuracy;
    sevenStepReturn(epsilonIndex) = metrics.sevenStepReturn;
    baselineTransportError(epsilonIndex) = ...
        metrics.baselineTransportError;
    lowModulusTransportError(epsilonIndex) = ...
        metrics.lowModulusTransportError;
    transportErrorEnrichment(epsilonIndex) = ...
        metrics.transportErrorEnrichment;
    metricData{epsilonIndex} = metrics;
end

spectralFigureIndex = find( ...
    abs(epsilonSweep - spectralFigureEpsilon) < ...
    10*eps(spectralFigureEpsilon), 1);
phaseFigureIndex = find( ...
    abs(epsilonSweep - phaseFigureEpsilon) < ...
    10*eps(phaseFigureEpsilon), 1);
assert(~isempty(spectralFigureIndex), ...
    'spectralFigureEpsilon must be a member of epsilonSweep.')
assert(~isempty(phaseFigureIndex), ...
    'phaseFigureEpsilon must be a member of epsilonSweep.')
phaseMetrics = metricData{phaseFigureIndex};
assert(phaseMetrics.regionShift == 3, ...
    'The oriented wave packet must advance three of the seven regions.')

robustnessTable = table( ...
    epsilonSweep(:), peakLocation/pi, targetMismatch, ...
    phaseConcentration, medianAbsoluteResidual, transitionAccuracy, ...
    sevenStepReturn, baselineTransportError, ...
    lowModulusTransportError, transportErrorEnrichment, ...
    'VariableNames', { ...
        'epsilon', 'thetaMaxOverPi', 'eta', 'rhoPhase', ...
        'medianAbsResidual', 'pPlus3', 'p7', 'baselineError', ...
        'lowModulusError', 'errorEnrichment'});

fprintf('\nSeven-region transport diagnostics\n')
fprintf('----------------------------------\n')
disp(robustnessTable)

%% Locate the three positive-frequency spectral peaks

companionTargets = [2*pi/7, 4*pi/7, 6*pi/7];
companionLocation = nan(numberOfEpsilonValues, 3);
companionHeight = nan(numberOfEpsilonValues, 3);

fprintf('Local maxima near 2*pi/7, 4*pi/7, and 6*pi/7:\n')
for epsilonIndex = 1:numberOfEpsilonValues
    for targetIndex = 1:numel(companionTargets)
        [companionLocation(epsilonIndex, targetIndex), ...
                companionHeight(epsilonIndex, targetIndex)] = ...
            local_maximum( ...
                thetaGrid, spectralDensitySweep(:, epsilonIndex), ...
                companionTargets(targetIndex), 0.06*pi);
    end
    fprintf('  epsilon = %.2f: theta/pi = [%7.4f %7.4f %7.4f]\n', ...
        epsilonSweep(epsilonIndex), ...
        companionLocation(epsilonIndex, :)/pi)
end

%% Check spectral robustness across scalar observables

observableNames = {'x', 'dx/dt', 'x+dx/dt'};
observableData = [x(:, 1), x(:, 2), x(:, 1) + x(:, 2)];
observableSpectra = nan(numel(thetaGrid), numel(observableNames));
observablePeakLocation = nan(1, numel(observableNames));
observablePeakHeight = nan(1, numel(observableNames));

if runObservableRobustness
    fprintf('\nObservable robustness at epsilon = %.2f\n', ...
        phaseFigureEpsilon)
    for observableIndex = 1:numel(observableNames)
        standardizedObservable = standardize_observable( ...
            observableData(:, observableIndex));
        [PsiXObservable, PsiYObservable] = build_delay_dictionary( ...
            standardizedObservable, delayDimension);
        [~, observableDensity] = run_rigged_dmd( ...
            PsiXObservable, PsiYObservable, uniformWeights, ...
            phaseFigureEpsilon, thetaTarget, kernelOrder, ...
            gCoefficients, thetaGrid);

        observableSpectra(:, observableIndex) = ...
            real(observableDensity(:));
        [observablePeakLocation(observableIndex), ...
                observablePeakHeight(observableIndex)] = ...
            local_maximum( ...
                thetaGrid, ...
                observableSpectra(:, observableIndex), ...
                thetaTarget, 0.06*pi);

        fprintf('  %-8s: peak theta/pi = %.6f, height = %.6e\n', ...
            observableNames{observableIndex}, ...
            observablePeakLocation(observableIndex)/pi, ...
            observablePeakHeight(observableIndex))
    end
end

%% Prepare the displayed seven-region partition

wavePacket = phaseMetrics.wavePacket;
phase = phaseMetrics.phase;
modulus = phaseMetrics.modulus;
logModulus = log(max(modulus, realmin));
region = phaseMetrics.region;
transitionCounts = phaseMetrics.transitionCounts;
transitionMatrix = phaseMetrics.transitionMatrix;

regionColours = lines(numberOfRegions);
regionLabels = arrayfun( ...
    @(regionIndex) sprintf('$R_%d$', regionIndex), ...
    1:numberOfRegions, 'UniformOutput', false);

%% Spectral-density figure at epsilon = 0.25

spectralFigure = figure( ...
    'Name', 'Duffing Regime II spectral density', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);

plot(thetaGrid/pi, spectralDensitySweep(:, spectralFigureIndex), ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 2.2)
hold on
for targetIndex = 1:numel(companionTargets)
    xline(companionTargets(targetIndex)/pi, 'k:', ...
        'LineWidth', 1.1, 'HandleVisibility', 'off');
end
xlabel('$\theta/\pi$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('Spectral density', ...
    'Interpreter', 'latex', 'FontSize', 20)
set(gca, 'FontSize', 16, ...
    'XTick', [0, 2/7, 4/7, 6/7, 1], ...
    'XTickLabel', {'$0$', '$2/7$', '$4/7$', '$6/7$', '$1$'}, ...
    'TickLabelInterpreter', 'latex')
xlim([0 1])
grid on
box on

%% Seven-region phase-partition figure

phaseFigure = figure( ...
    'Name', 'Duffing Regime II phase partition', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);

scatter(sectionData(:, 1), sectionData(:, 2), ...
    12, region, 'filled')
xlabel('$x$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\dot{x}$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
colormap(phaseFigure, regionColours)
caxis([0.5, numberOfRegions + 0.5])
phaseColourBar = colorbar;
phaseColourBar.Ticks = 1:numberOfRegions;
phaseColourBar.TickLabels = regionLabels;
phaseColourBar.TickLabelInterpreter = 'latex';
axis tight
grid on
box on

%% One-step transport figure

transportFigure = figure( ...
    'Name', 'Duffing Regime II one-step transport', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);

sourceRegion = region(1:end-1);
nextState = sectionData(2:end, :);
xLimits = [min(sectionData(:, 1)), max(sectionData(:, 1))];
velocityLimits = [min(sectionData(:, 2)), max(sectionData(:, 2))];

scatter(nextState(:, 1), nextState(:, 2), ...
    12, sourceRegion, 'filled')
xlabel('$x$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\dot{x}$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
xlim(xLimits)
ylim(velocityLimits)
caxis([0.5, numberOfRegions + 0.5])
colormap(transportFigure, regionColours)
transportColourBar = colorbar;
transportColourBar.Ticks = 1:numberOfRegions;
transportColourBar.TickLabels = regionLabels;
transportColourBar.TickLabelInterpreter = 'latex';
grid on
box on

%% Archive the analysis and write its metadata report

analysis = struct();
analysis.generatedOn = datestr(now, 30);
analysis.observable = 'standardized position x';
analysis.observableMean = positionMean;
analysis.observableScale = positionScale;
analysis.N = N;
analysis.delayDimension = delayDimension;
analysis.nSnapshotPairs = nSnapshotPairs;
analysis.uniformWeight = uniformWeights(1);
analysis.uniformWeightSum = uniformWeightSum;
analysis.epsilonSweep = epsilonSweep;
analysis.spectralFigureEpsilon = spectralFigureEpsilon;
analysis.phaseFigureEpsilon = phaseFigureEpsilon;
analysis.kernelOrder = kernelOrder;
analysis.thetaTarget = thetaTarget;
analysis.numberOfRegions = numberOfRegions;
analysis.reliabilityQuantile = reliabilityQuantile;
analysis.lowModulusQuantile = lowModulusQuantile;
analysis.robustnessTable = robustnessTable;
analysis.peakLocation = peakLocation;
analysis.peakHeight = peakHeight;
analysis.meanPhaseAdvance = meanPhaseAdvance;
analysis.targetMismatch = targetMismatch;
analysis.phaseConcentration = phaseConcentration;
analysis.medianAbsoluteResidual = medianAbsoluteResidual;
analysis.transitionAccuracy = transitionAccuracy;
analysis.sevenStepReturn = sevenStepReturn;
analysis.baselineTransportError = baselineTransportError;
analysis.lowModulusTransportError = lowModulusTransportError;
analysis.transportErrorEnrichment = transportErrorEnrichment;
analysis.companionTargets = companionTargets;
analysis.companionLocation = companionLocation;
analysis.companionHeight = companionHeight;
analysis.observableNames = observableNames;
analysis.observablePeakLocation = observablePeakLocation;
analysis.observablePeakHeight = observablePeakHeight;
analysis.transitionCounts = transitionCounts;
analysis.transitionMatrix = transitionMatrix;
analysis.riggedDMDImplementation = which('riggedDMD');

save(cacheFile, ...
    'x', 'sectionData', 'generationMetadata', 'analysis', ...
    'thetaGrid', 'spectralDensitySweep', ...
    'epsilonSweep', 'spectralFigureEpsilon', 'phaseFigureEpsilon', ...
    'metricData', 'phaseMetrics', ...
    'wavePacket', 'phase', 'modulus', 'logModulus', 'region', ...
    'transitionCounts', 'transitionMatrix', ...
    'companionLocation', 'companionHeight', ...
    'observableNames', 'observableSpectra', ...
    'observablePeakLocation', 'observablePeakHeight', '-v7')

write_metadata_report(metadataFile, generationMetadata, analysis)

if saveFigures
    save_figure_pdf(spectralFigure, fullfile( ...
        outputDirectory, 'duffing_regime2_spectral_density.pdf'))
    save_figure_pdf(phaseFigure, fullfile( ...
        outputDirectory, 'duffing_regime2_phase_partition.pdf'))
    save_figure_pdf(transportFigure, fullfile( ...
        outputDirectory, 'duffing_regime2_one_step_transport.pdf'))
end

fprintf('\nSaved the section data and Regime II analysis to %s\n', cacheFile)
fprintf('Saved provenance and diagnostics to %s\n', metadataFile)
fprintf('Finished writing the three Regime II PDF figures.\n')

%% Local functions

function [x, metadata, generatedNow] = ...
        load_or_generate_section_data(cacheFile, config, forceRegenerate)
    %LOAD_OR_GENERATE_SECTION_DATA Use only a documented matching cache.

    generatedNow = false;
    cacheIsUsable = false;
    reason = 'cache file not found';

    if isfile(cacheFile) && ~forceRegenerate
        loaded = load(cacheFile);
        [cacheIsUsable, reason] = validate_section_cache(loaded, config);
        if cacheIsUsable
            x = loaded.x;
            metadata = loaded.generationMetadata;
            return
        end
    elseif forceRegenerate
        reason = 'regeneration explicitly requested';
    end

    fprintf('Generating the event-located Duffing section (%s).\n', reason)
    [x, metadata] = generate_section_data(config);
    generationMetadata = metadata; 
    save(cacheFile, 'x', 'generationMetadata', '-v7')
    generatedNow = true;
end

function [valid, reason] = validate_section_cache(loaded, config)
    %VALIDATE_SECTION_CACHE Check data shape, provenance, and parameters.

    valid = false;
    reason = '';

    if ~isfield(loaded, 'x')
        reason = 'cache does not contain x';
        return
    end
    if ~isnumeric(loaded.x) || size(loaded.x, 2) ~= 2 || ...
            any(~isfinite(loaded.x(:)))
        reason = 'cache variable x is not a finite N-by-2 array';
        return
    end
    if ~isfield(loaded, 'generationMetadata')
        reason = 'legacy cache does not contain generation metadata';
        return
    end

    metadata = loaded.generationMetadata;
    requiredFields = { ...
        'schemaVersion', 'parameters', 'initialCondition', ...
        'integrationInterval', 'outputStep', ...
        'numberOfOriginalOutputTimes', 'solver', 'relativeTolerance', ...
        'absoluteTolerance', 'maximumStep', 'refine', ...
        'integrationBlockLength', 'discardedCrossings', ...
        'forcingPhase', 'phasePeriod', 'forcingPeriod', ...
        'eventDirection', 'crossingRule', 'eventLocationRule', ...
        'retainedEventTimes', 'maximumEventPhaseResidual', 'N'};

    for fieldIndex = 1:numel(requiredFields)
        if ~isfield(metadata, requiredFields{fieldIndex})
            reason = ['cache lacks provenance field ' ...
                requiredFields{fieldIndex}];
            return
        end
    end

    parametersMatch = ...
        metadata.parameters.alpha == config.alpha && ...
        metadata.parameters.beta == config.beta && ...
        metadata.parameters.delta == config.delta && ...
        metadata.parameters.omega == config.omega && ...
        metadata.parameters.gamma == config.gamma;

    settingsMatch = ...
        metadata.schemaVersion == config.schemaVersion && ...
        isequal(metadata.initialCondition(:), ...
            config.initialCondition(:)) && ...
        isequal(metadata.integrationInterval, ...
            config.integrationInterval) && ...
        metadata.outputStep == config.outputStep && ...
        metadata.numberOfOriginalOutputTimes == ...
            config.numberOfOriginalOutputTimes && ...
        strcmp(metadata.solver, config.solver) && ...
        metadata.relativeTolerance == config.relativeTolerance && ...
        isequal(metadata.absoluteTolerance, ...
            config.absoluteTolerance) && ...
        metadata.maximumStep == config.maximumStep && ...
        metadata.refine == config.refine && ...
        metadata.integrationBlockLength == ...
            config.integrationBlockLength && ...
        metadata.discardedCrossings == config.discardedCrossings && ...
        metadata.forcingPhase == config.forcingPhase && ...
        metadata.phasePeriod == config.phasePeriod && ...
        metadata.forcingPeriod == config.forcingPeriod && ...
        metadata.eventDirection == config.eventDirection;

    if ~(parametersMatch && settingsMatch)
        reason = 'cache metadata do not match the requested settings';
        return
    end
    if metadata.N ~= size(loaded.x, 1) || ...
            numel(metadata.retainedEventTimes) ~= metadata.N
        reason = 'cache metadata and section-array lengths disagree';
        return
    end

    valid = true;
    reason = 'documented cache matches the requested settings';
end

function [x, metadata] = generate_section_data(config)
    %GENERATE_SECTION_DATA Integrate in blocks and retain event roots only.

    eventFunction = @(time, state) ...
        stroboscopic_event( ...
            time, state, config.forcingPhase, config.eventDirection);
    odeOptions = odeset( ...
        'RelTol', config.relativeTolerance, ...
        'AbsTol', config.absoluteTolerance, ...
        'MaxStep', config.maximumStep, ...
        'Refine', config.refine, ...
        'Events', eventFunction);

    tStart = config.integrationInterval(1);
    tFinal = config.integrationInterval(2);
    state = config.initialCondition(:);
    blockStart = tStart;
    blockIndex = 0;
    numberOfBlocks = ceil((tFinal - tStart)/ ...
        config.integrationBlockLength);

    eventTimeCells = cell(numberOfBlocks, 1);
    eventStateCells = cell(numberOfBlocks, 1);

    fprintf(['Integrating on [%.16g, %.16g] with %s in %d blocks.\n'], ...
        tStart, tFinal, config.solver, numberOfBlocks)

    while blockStart < tFinal
        blockIndex = blockIndex + 1;
        blockEnd = min( ...
            blockStart + config.integrationBlockLength, tFinal);

        [~, blockSolution, eventTimes, eventStates] = ode45( ...
            @(time, currentState) duffing_rhs( ...
                time, currentState, config), ...
            [blockStart, blockEnd], state, odeOptions);

        state = blockSolution(end, :).';
        eventTimeCells{blockIndex} = eventTimes;
        eventStateCells{blockIndex} = eventStates;
        blockStart = blockEnd;

        if mod(blockIndex, 20) == 0 || blockIndex == numberOfBlocks
            fprintf('  completed block %d of %d (t = %.16g)\n', ...
                blockIndex, numberOfBlocks, blockEnd)
        end
    end

    rawEventTimes = vertcat(eventTimeCells{:});
    rawEventStates = vertcat(eventStateCells{:});

    assert(~isempty(rawEventTimes), ...
        'No stroboscopic events were located on the integration interval.')

    [rawEventTimes, sortOrder] = sort(rawEventTimes);
    rawEventStates = rawEventStates(sortOrder, :);

    timeTolerance = max(1e-10, 100*eps(max(abs(rawEventTimes))));
    afterInitialPoint = rawEventTimes > tStart + timeTolerance;
    rawEventTimes = rawEventTimes(afterInitialPoint);
    rawEventStates = rawEventStates(afterInitialPoint, :);

    uniqueEvent = [true; diff(rawEventTimes) > timeTolerance];
    rawEventTimes = rawEventTimes(uniqueEvent);
    rawEventStates = rawEventStates(uniqueEvent, :);

    phaseResidual = atan2( ...
        sin(rawEventStates(:, 3) - config.forcingPhase), ...
        cos(rawEventStates(:, 3) - config.forcingPhase));

    assert(config.omega > 0, ...
        'The positive-direction event rule requires omega>0.')
    assert(numel(rawEventTimes) > config.discardedCrossings, ...
        'All stroboscopic crossings were removed by the transient filter.')

    retainedIndices = (config.discardedCrossings + 1): ...
        numel(rawEventTimes);
    retainedEventTimes = rawEventTimes(retainedIndices);
    retainedEventStates = rawEventStates(retainedIndices, :);
    x = retainedEventStates(:, 1:2);

    initialRelativePhase = ...
        config.initialCondition(3) - config.forcingPhase;
    finalRelativePhase = initialRelativePhase + ...
        config.omega*(tFinal - tStart);
    expectedRawCrossings = ...
        floor(finalRelativePhase/config.phasePeriod) - ...
        floor(initialRelativePhase/config.phasePeriod);
    assert(numel(rawEventTimes) == expectedRawCrossings, ...
        ['Expected %d positive-phase events but located %d. ' ...
         'Check the event rule and integration settings.'], ...
        expectedRawCrossings, numel(rawEventTimes))

    metadata = struct();
    metadata.schemaVersion = config.schemaVersion;
    metadata.createdOn = datestr(now, 30);
    metadata.system = config.system;
    metadata.parameters = struct( ...
        'alpha', config.alpha, ...
        'beta', config.beta, ...
        'delta', config.delta, ...
        'omega', config.omega, ...
        'gamma', config.gamma);
    metadata.initialCondition = config.initialCondition;
    metadata.integrationInterval = config.integrationInterval;
    metadata.outputStep = config.outputStep;
    metadata.numberOfOriginalOutputTimes = ...
        config.numberOfOriginalOutputTimes;
    metadata.solver = config.solver;
    metadata.relativeTolerance = config.relativeTolerance;
    metadata.absoluteTolerance = config.absoluteTolerance;
    metadata.maximumStep = config.maximumStep;
    metadata.refine = config.refine;
    metadata.integrationBlockLength = config.integrationBlockLength;
    metadata.forcingPhase = config.forcingPhase;
    metadata.phasePeriod = config.phasePeriod;
    metadata.forcingPeriod = config.forcingPeriod;
    metadata.eventDirection = config.eventDirection;
    metadata.crossingRule = config.crossingRule;
    metadata.eventLocationRule = config.eventLocationRule;
    metadata.rawEventCount = numel(rawEventTimes);
    metadata.discardedCrossings = config.discardedCrossings;
    if config.discardedCrossings == 0
        metadata.discardedTransientInterval = [tStart, tStart];
    else
        metadata.discardedTransientInterval = [ ...
            tStart, retainedEventTimes(1)];
    end
    metadata.retainedEventTimes = retainedEventTimes;
    metadata.retainedEventStates = retainedEventStates;
    metadata.maximumEventPhaseResidual = max(abs(phaseResidual));
    metadata.N = size(x, 1);
end

function [value, isTerminal, direction] = ...
        stroboscopic_event(~, state, forcingPhase, eventDirection)
    %STROBOSCOPIC_EVENT Locate phi=0 mod 2*pi in positive direction.

    value = sin(state(3) - forcingPhase);
    isTerminal = 0;
    direction = eventDirection;
end

function derivative = duffing_rhs(~, state, config)
    %DUFFING_RHS Autonomous phase formulation of the forced Duffing ODE.

    derivative = [ ...
        state(2); ...
        -config.delta*state(2) ...
        - config.alpha*state(1) ...
        - config.beta*state(1)^3 ...
        + config.gamma*cos(state(3)); ...
        config.omega];
end

function [standardized, center, scale] = ...
        standardize_observable(observable)
    %STANDARDIZE_OBSERVABLE Center and scale a scalar trajectory.

    observable = observable(:);
    assert(all(isfinite(observable)), ...
        'The observable contains NaN or Inf values.')

    center = mean(observable);
    scale = std(observable);
    assert(scale > 0, 'The observable has zero variance.')
    standardized = (observable - center)/scale;
end

function [PsiX, PsiY] = build_delay_dictionary(observable, dimension)
    %BUILD_DELAY_DICTIONARY Construct paired delay-coordinate matrices.

    observable = observable(:);
    nRows = numel(observable) - dimension;
    assert(nRows > 0, ...
        'The delay dimension must be smaller than the data length.')

    rowOffsets = (1:nRows).';
    delayOffsets = 0:(dimension - 1);
    indices = rowOffsets + delayOffsets;

    PsiX = observable(indices);
    PsiY = observable(indices + 1);
end

function [mode, spectralDensity] = run_rigged_dmd( ...
        PsiX, PsiY, weight, epsilon, theta, kernelOrder, ...
        gCoefficients, thetaGrid)
    %RUN_RIGGED_DMD Evaluate one wave packet and its spectral density.

    [mode, spectralDensity] = riggedDMD( ...
        PsiX, PsiY, weight, epsilon, theta, [], ...
        'order', kernelOrder, ...
        'g_coeffs', gCoefficients, ...
        'TH2', thetaGrid);
    mode = squeeze(mode);
end

function metrics = transport_metrics( ...
        rawWavePacket, thetaTarget, reliabilityQuantile, ...
        lowModulusQuantile, numberOfRegions)
    %TRANSPORT_METRICS Compute the Regime II phase-transport diagnostics.

    rawWavePacket = rawWavePacket(:);
    modulus = abs(rawWavePacket);
    modulusThreshold = quantile(modulus, reliabilityQuantile);
    reliableStep = ...
        modulus(1:end-1) > modulusThreshold & ...
        modulus(2:end) > modulusThreshold;
    assert(any(reliableStep), ...
        'The reliability filter removed every one-step pair.')

    rawPhase = angle(rawWavePacket);
    rawIncrement = wrap_to_pi( ...
        rawPhase(2:end) - rawPhase(1:end-1));
    rawMeanAdvance = circular_mean(rawIncrement(reliableStep));

    positiveMismatch = abs(wrap_to_pi(rawMeanAdvance - thetaTarget));
    negativeMismatch = abs(wrap_to_pi(rawMeanAdvance + thetaTarget));
    conjugated = negativeMismatch < positiveMismatch;

    if conjugated
        wavePacket = conj(rawWavePacket);
    else
        wavePacket = rawWavePacket;
    end

    phase = angle(wavePacket);
    phaseIncrement = wrap_to_pi(phase(2:end) - phase(1:end-1));
    meanPhaseAdvance = circular_mean(phaseIncrement(reliableStep));
    phaseResidual = wrap_to_pi(phaseIncrement - thetaTarget);
    reliableResidual = phaseResidual(reliableStep);

    targetMismatch = abs(wrap_to_pi(meanPhaseAdvance - thetaTarget));
    phaseConcentration = abs(mean(exp(1i*reliableResidual)));
    medianAbsoluteResidual = median(abs(reliableResidual));

    regionWidth = 2*pi/numberOfRegions;
    region = floor((phase + pi)/regionWidth) + 1;
    region = min(max(region, 1), numberOfRegions);
    regionShift = round(thetaTarget/regionWidth);

    sourceRegion = region(1:end-1);
    observedDestination = region(2:end);
    predictedDestination = mod( ...
        sourceRegion - 1 + regionShift, numberOfRegions) + 1;

    [transitionCounts, transitionMatrix] = ...
        empirical_transition_matrix(region, numberOfRegions);

    followsTransportRule = ...
        observedDestination == predictedDestination;
    transitionAccuracy = mean(followsTransportRule(reliableStep));

    sevenStepReturn = mean( ...
        region(1:end-numberOfRegions) == ...
        region(1 + numberOfRegions:end));

    stepModulus = min( ...
        modulus(1:end-1), modulus(2:end));
    retainedStepModulus = stepModulus(reliableStep);
    lowModulusThreshold = quantile( ...
        retainedStepModulus, lowModulusQuantile);
    lowModulusStep = ...
        reliableStep & stepModulus <= lowModulusThreshold;
    assert(any(lowModulusStep), ...
        'The retained low-modulus set is empty.')

    transportError = ~followsTransportRule;
    baselineTransportError = mean(transportError(reliableStep));
    lowModulusTransportError = mean(transportError(lowModulusStep));
    assert(baselineTransportError > 0, ...
        'Transport-error enrichment is undefined when no errors occur.')
    transportErrorEnrichment = ...
        lowModulusTransportError/baselineTransportError;

    metrics = struct();
    metrics.wavePacket = wavePacket;
    metrics.conjugated = conjugated;
    metrics.phase = phase;
    metrics.modulus = modulus;
    metrics.modulusThreshold = modulusThreshold;
    metrics.phaseIncrement = phaseIncrement;
    metrics.phaseResidual = phaseResidual;
    metrics.reliableStep = reliableStep;
    metrics.reliableFraction = mean(reliableStep);
    metrics.meanPhaseAdvance = meanPhaseAdvance;
    metrics.targetMismatch = targetMismatch;
    metrics.phaseConcentration = phaseConcentration;
    metrics.medianAbsoluteResidual = medianAbsoluteResidual;
    metrics.region = region;
    metrics.regionShift = regionShift;
    metrics.transitionCounts = transitionCounts;
    metrics.transitionMatrix = transitionMatrix;
    metrics.transitionAccuracy = transitionAccuracy;
    metrics.sevenStepReturn = sevenStepReturn;
    metrics.stepModulus = stepModulus;
    metrics.lowModulusThreshold = lowModulusThreshold;
    metrics.lowModulusStep = lowModulusStep;
    metrics.transportError = transportError;
    metrics.baselineTransportError = baselineTransportError;
    metrics.lowModulusTransportError = lowModulusTransportError;
    metrics.transportErrorEnrichment = transportErrorEnrichment;
end

function meanAngle = circular_mean(angleValues)
    %CIRCULAR_MEAN Return the argument of the mean unit phasor.

    assert(~isempty(angleValues), ...
        'Cannot compute a circular mean of an empty array.')
    meanAngle = angle(mean(exp(1i*angleValues)));
end

function wrapped = wrap_to_pi(angleValues)
    %WRAP_TO_PI Map angles to [-pi,pi).

    wrapped = mod(angleValues + pi, 2*pi) - pi;
end

function [counts, probabilities] = ...
        empirical_transition_matrix(region, numberOfRegions)
    %EMPIRICAL_TRANSITION_MATRIX Count all available one-step pairs.

    source = region(1:end-1);
    destination = region(2:end);
    counts = accumarray( ...
        [source, destination], 1, ...
        [numberOfRegions, numberOfRegions]);

    rowTotals = sum(counts, 2);
    probabilities = zeros(size(counts));
    occupiedRows = rowTotals > 0;
    probabilities(occupiedRows, :) = ...
        counts(occupiedRows, :)./rowTotals(occupiedRows);
end

function [location, height] = ...
        local_maximum(thetaGrid, density, target, halfWidth)
    %LOCAL_MAXIMUM Find the sampled density maximum near a target angle.

    thetaGrid = thetaGrid(:);
    density = density(:);
    localMask = ...
        thetaGrid >= target - halfWidth & ...
        thetaGrid <= target + halfWidth & ...
        isfinite(density);

    assert(any(localMask), ...
        'No finite spectral-density value lies in the search window.')
    localIndices = find(localMask);
    [height, relativeIndex] = max(density(localMask));
    location = thetaGrid(localIndices(relativeIndex));
end

function write_metadata_report(fileName, metadata, analysis)
    %WRITE_METADATA_REPORT Write provenance and diagnostics as plain text.

    fileID = fopen(fileName, 'w');
    assert(fileID >= 0, 'Could not open %s for writing.', fileName)
    closeFile = onCleanup(@() fclose(fileID)); 

    fprintf(fileID, 'Forced Duffing Regime II stroboscopic archive\n');
    fprintf(fileID, 'Created: %s\n', metadata.createdOn);
    fprintf(fileID, 'Analyzed: %s\n\n', analysis.generatedOn);

    fprintf(fileID, ...
        'Parameters: alpha=%.16g, beta=%.16g, delta=%.16g, ', ...
        metadata.parameters.alpha, metadata.parameters.beta, ...
        metadata.parameters.delta);
    fprintf(fileID, 'omega=%.16g, gamma=%.16g\n', ...
        metadata.parameters.omega, metadata.parameters.gamma);
    fprintf(fileID, ...
        'Initial condition [x,v,phi]: [%.16g, %.16g, %.16g]\n', ...
        metadata.initialCondition);
    fprintf(fileID, 'Integration interval: [%.16g, %.16g]\n', ...
        metadata.integrationInterval);
    fprintf(fileID, 'Original requested output spacing: %.16g\n', ...
        metadata.outputStep);
    fprintf(fileID, 'Original requested output count: %.16g\n', ...
        metadata.numberOfOriginalOutputTimes);
    fprintf(fileID, 'Solver: %s\n', metadata.solver);
    fprintf(fileID, 'Relative tolerance: %.16g\n', ...
        metadata.relativeTolerance);
    fprintf(fileID, 'Absolute tolerances: [%.16g, %.16g, %.16g]\n', ...
        metadata.absoluteTolerance);
    fprintf(fileID, 'Maximum solver step: %.16g\n', ...
        metadata.maximumStep);
    fprintf(fileID, 'Solver output refinement: %d\n', metadata.refine);
    fprintf(fileID, 'Integration block length: %.16g\n', ...
        metadata.integrationBlockLength);
    fprintf(fileID, 'Forcing phase at section: %.16g rad\n', ...
        metadata.forcingPhase);
    fprintf(fileID, 'Forcing period: %.16g\n', metadata.forcingPeriod);
    fprintf(fileID, 'Crossing rule: %s\n', metadata.crossingRule);
    fprintf(fileID, 'Event location: %s\n', metadata.eventLocationRule);
    fprintf(fileID, 'Raw event count: %d\n', metadata.rawEventCount);
    fprintf(fileID, 'Discarded initial crossings: %d\n', ...
        metadata.discardedCrossings);
    fprintf(fileID, 'Discarded transient interval: [%.16g, %.16g)\n', ...
        metadata.discardedTransientInterval);
    fprintf(fileID, 'Retained event interval: [%.16g, %.16g]\n', ...
        metadata.retainedEventTimes(1), metadata.retainedEventTimes(end));
    fprintf(fileID, 'Maximum event phase residual: %.16e rad\n\n', ...
        metadata.maximumEventPhaseResidual);

    fprintf(fileID, 'Retained section points N: %d\n', analysis.N);
    fprintf(fileID, 'Observable: %s\n', analysis.observable);
    fprintf(fileID, 'Observable mean: %.16g\n', analysis.observableMean);
    fprintf(fileID, 'Observable standard deviation: %.16g\n', ...
        analysis.observableScale);
    fprintf(fileID, 'Delay dimension d: %d\n', analysis.delayDimension);
    fprintf(fileID, 'Hankel snapshot pairs N-d: %d\n', ...
        analysis.nSnapshotPairs);
    fprintf(fileID, 'Uniform quadrature weight: %.16g\n', ...
        analysis.uniformWeight);
    fprintf(fileID, 'Sum of weights: %.16g\n', ...
        analysis.uniformWeightSum);
    fprintf(fileID, 'Target theta: %.16g = %.16g*pi\n', ...
        analysis.thetaTarget, analysis.thetaTarget/pi);
    fprintf(fileID, 'Reliability quantile: %.16g\n', ...
        analysis.reliabilityQuantile);
    fprintf(fileID, 'Retained low-modulus quantile: %.16g\n\n', ...
        analysis.lowModulusQuantile);

    fprintf(fileID, ['epsilon theta_max/pi eta rho_phase med_abs_r ' ...
        'p_plus3 p7 baseline_error low_modulus_error E_err\n']);
    diagnostics = analysis.robustnessTable;
    for rowIndex = 1:height(diagnostics)
        fprintf(fileID, ...
            ['%7.4f %12.8f %12.8e %12.8f %12.8f ' ...
             '%12.8f %12.8f %12.8f %12.8f %12.8f\n'], ...
            diagnostics.epsilon(rowIndex), ...
            diagnostics.thetaMaxOverPi(rowIndex), ...
            diagnostics.eta(rowIndex), ...
            diagnostics.rhoPhase(rowIndex), ...
            diagnostics.medianAbsResidual(rowIndex), ...
            diagnostics.pPlus3(rowIndex), ...
            diagnostics.p7(rowIndex), ...
            diagnostics.baselineError(rowIndex), ...
            diagnostics.lowModulusError(rowIndex), ...
            diagnostics.errorEnrichment(rowIndex));
    end

    fprintf(fileID, '\nOriented circular mean phase advance by epsilon:\n');
    for rowIndex = 1:height(diagnostics)
        fprintf(fileID, '  epsilon=%.4f mean_advance/pi=%.10f\n', ...
            diagnostics.epsilon(rowIndex), ...
            analysis.meanPhaseAdvance(rowIndex)/pi);
    end

    spectralFigureIndex = find( ...
        abs(diagnostics.epsilon - analysis.spectralFigureEpsilon) < ...
        10*eps(analysis.spectralFigureEpsilon), 1);
    fprintf(fileID, '\nSpectral-figure epsilon: %.16g\n', ...
        analysis.spectralFigureEpsilon);
    fprintf(fileID, 'Phase/transport-figure epsilon: %.16g\n', ...
        analysis.phaseFigureEpsilon);
    fprintf(fileID, ...
        'Spectral-figure peak locations theta/pi: %.8f %.8f %.8f\n', ...
        analysis.companionLocation(spectralFigureIndex, :)/pi);
    fprintf(fileID, ...
        'Spectral-figure peak heights: %.8e %.8e %.8e\n', ...
        analysis.companionHeight(spectralFigureIndex, :));

    if all(isfinite(analysis.observablePeakLocation))
        fprintf(fileID, '\nObservable peaks at phase-figure epsilon:\n');
        for observableIndex = 1:numel(analysis.observableNames)
            fprintf(fileID, '  %-8s theta/pi=%.8f height=%.8e\n', ...
                analysis.observableNames{observableIndex}, ...
                analysis.observablePeakLocation(observableIndex)/pi, ...
                analysis.observablePeakHeight(observableIndex));
        end
    end

    fprintf(fileID, ...
        '\nDisplayed transition counts (rows source, columns destination):\n');
    write_matrix(fileID, analysis.transitionCounts, '%10d');
    fprintf(fileID, 'Displayed row-normalized transition matrix:\n');
    write_matrix(fileID, analysis.transitionMatrix, '%12.8f');
    fprintf(fileID, '\nriggedDMD implementation: %s\n', ...
        analysis.riggedDMDImplementation);
end

function write_matrix(fileID, matrix, format)
    %WRITE_MATRIX Write a small numeric matrix with aligned columns.

    for rowIndex = 1:size(matrix, 1)
        for columnIndex = 1:size(matrix, 2)
            fprintf(fileID, format, matrix(rowIndex, columnIndex));
        end
        fprintf(fileID, '\n');
    end
end

function save_figure_pdf(figHandle, fileName)
    %SAVE_FIGURE_PDF Save a figure to a tightly sized vector PDF.

    [fileDirectory, ~, extension] = fileparts(fileName);
    if isempty(extension)
        fileName = [fileName '.pdf'];
    elseif ~strcmpi(extension, '.pdf')
        error('Figure output must have a .pdf extension.')
    end

    assert(isempty(fileDirectory) || isfolder(fileDirectory), ...
        'The figure output directory does not exist: %s', fileDirectory)

    set(figHandle, 'PaperUnits', 'centimeters')
    set(figHandle, 'Units', 'centimeters')
    position = get(figHandle, 'Position');
    set(figHandle, 'PaperSize', [position(3), position(4)])
    set(figHandle, 'PaperPositionMode', 'manual')
    set(figHandle, 'PaperPosition', [0, 0, position(3), position(4)])
    print(figHandle, fileName, '-dpdf', '-painters')
end
