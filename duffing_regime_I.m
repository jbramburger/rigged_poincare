% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE FORCED DUFFING POINCARE MAP: REGIME I
% -------------------------------------------------------------------------
%
% This script generates a documented, event-located stroboscopic section
% for the forced Duffing oscillator and applies rigged Dynamic Mode
% Decomposition (riggedDMD) to its dominant three-fold organization.
%
% The autonomous phase formulation is
%
%   xdot     = v,
%   vdot     = -delta*v - alpha*x - beta*x^3 + gamma*cos(phi),
%   phidot   = omega.
%
% Stroboscopic intersections are roots of sin(phi)=0 with positive
% direction, equivalently phi=0 mod 2*pi. MATLAB's continuous event
% interpolation is used to locate each crossing. The trajectory is
% integrated in blocks to bound memory use; the final state of each block
% is passed unchanged to the next block.
%
% The analysis uses the standardized position observable, d=20 delays,
% epsilon=0.25, and theta=2*pi/3. It computes the spectral density, finite
% regularized wave packet, three phase regions, empirical transition
% matrix, phase-advance diagnostics, and an observable-robustness check.
% The archived analysis is read by duffing_utils/duffing_upo_plot.m, which
% produces the six phase/modulus panels in the paper without repeating the
% UPO shooting calculation.
%
% REQUIRED CODE
%   main_routines/riggedDMD.m and its dependencies
%
% OUTPUT DATA
%   duffing_section_I_cache.mat
%       Event-located section data, provenance, wave-packet values,
%       phase regions, spectra, and diagnostics
%   duffing_section_I_metadata.txt
%       Human-readable generation and analysis metadata
%
% All output files are written beside this script. No output directory is
% created.
%
% OUTPUT FIGURES
%   duffing_regime1_spectral_density.pdf
%   duffing_regime1_observable_robustness.pdf
%   duffing_regime1_phase_partition.pdf
%   duffing_regime1_modulus.pdf
%
% PAPER PANELS
%   After this script finishes, run
%
%       run('duffing_utils/duffing_upo_plot.m')
%
%   to regenerate duffing_phase_1.pdf,...,duffing_phase_3.pdf and
%   duffing_mod_1.pdf,...,duffing_mod_3.pdf. Existing UPO data are reused.
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
cacheFile = fullfile(outputDirectory, 'duffing_section_I_cache.mat');
metadataFile = fullfile( ...
    outputDirectory, 'duffing_section_I_metadata.txt');
routinesDirectory = fullfile(scriptDirectory, 'main_routines');

forceRegenerateSectionData = false;
saveFigures = true;
runObservableRobustness = true;

%% Duffing model and section-generation parameters

config = struct();
config.schemaVersion = 2;
config.system = 'Forced Duffing oscillator, Regime I';
config.alpha = -1;
config.beta = 1;
config.delta = 0.3;
config.omega = 1.2;
config.gamma = 0.5;
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
config.discardedCrossings = 98;
config.forcingPhase = 0;
config.eventDirection = +1;
config.phasePeriod = 2*pi;
config.forcingPeriod = 2*pi/config.omega;
config.crossingRule = [ ...
    'Event-located root of sin(phi)=0 with positive direction, ' ...
    'excluding the initial point; equivalently phi=0 mod 2*pi.'];
config.eventLocationRule = [ ...
    'Root located by the continuous event interpolation used by ode45.'];

%% riggedDMD parameters

delayDimension = 20;
epsilon = 0.25;
kernelOrder = 2;
thetaTarget = 2*pi/3;
thetaGrid = linspace(-pi, pi, 1601);
reliabilityQuantile = 0.05;

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
expectedRetainedPoints = expectedRawCrossings - config.discardedCrossings;

assert(N == expectedRetainedPoints, ...
    'Expected %d retained samples, but the cache contains %d.', ...
    expectedRetainedPoints, N)
assert(N == 19000, ...
    'The documented Regime I computation must retain N=19000 samples.')

fprintf('\nDuffing Regime I section data\n')
fprintf('-----------------------------\n')
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

uniformWeights = ones(nSnapshotPairs, 1)/nSnapshotPairs;
assert(abs(sum(uniformWeights) - 1) < 1e-12, ...
    'The uniform quadrature weights do not sum to one.')

gCoefficients = zeros(delayDimension, 1);
gCoefficients(1) = 1;

fprintf('\nriggedDMD data\n')
fprintf('---------------\n')
fprintf('Observable = standardized position x\n')
fprintf('Delay dimension d = %d\n', delayDimension)
fprintf('Hankel snapshot pairs N-d = %d\n', nSnapshotPairs)
fprintf('Uniform quadrature weight = 1/%d = %.16g\n', ...
    nSnapshotPairs, uniformWeights(1))
fprintf('Sum of weights = %.16g\n', sum(uniformWeights))

%% Spectral density and finite wave packet

fprintf('Running riggedDMD with epsilon=%.4g and theta=%.10f.\n', ...
    epsilon, thetaTarget)

[wavePacketCoefficients, spectralDensity] = riggedDMD( ...
    PsiX, PsiY, uniformWeights, epsilon, thetaTarget, [], ...
    'order', kernelOrder, ...
    'g_coeffs', gCoefficients, ...
    'TH2', thetaGrid);

wavePacketCoefficients = squeeze(wavePacketCoefficients);
spectralDensity = real(squeeze(spectralDensity));
spectralDensity = spectralDensity(:).';

phi = PsiX*wavePacketCoefficients;
phase = angle(phi);
modulus = abs(phi);
logModulus = log(max(modulus, realmin));

[spectralPeakTheta, spectralPeakValue] = local_maximum( ...
    thetaGrid, spectralDensity, thetaTarget, 0.08*pi);

fprintf('Spectral peak near 2*pi/3: theta = %.10f = %.10f*pi\n', ...
    spectralPeakTheta, spectralPeakTheta/pi)
fprintf('Spectral peak value = %.10g\n', spectralPeakValue)

%% Phase advance and three-region partition

phaseIncrement = wrap_to_pi(phase(2:end) - phase(1:end-1));
modulusThreshold = quantile(modulus, reliabilityQuantile);
reliableStep = ...
    modulus(1:end-1) > modulusThreshold & ...
    modulus(2:end) > modulusThreshold;

assert(any(reliableStep), ...
    'No adjacent samples survived the low-modulus reliability filter.')

meanPhaseAdvance = angle(mean(exp(1i*phaseIncrement(reliableStep))));
positiveMismatch = abs(wrap_to_pi(meanPhaseAdvance - thetaTarget));
negativeMismatch = abs(wrap_to_pi(meanPhaseAdvance + thetaTarget));

if positiveMismatch <= negativeMismatch
    phaseSign = 1;
    phaseResidual = wrap_to_pi(phaseIncrement - thetaTarget);
    targetMismatch = positiveMismatch;
else
    phaseSign = -1;
    phaseResidual = wrap_to_pi(phaseIncrement + thetaTarget);
    targetMismatch = negativeMismatch;
end

phaseConcentration = abs(mean(exp(1i*phaseResidual(reliableStep))));

numberOfRegions = 3;
regionWidth = 2*pi/numberOfRegions;
region = floor((phase + pi)/regionWidth) + 1;
region = min(max(region, 1), numberOfRegions);

[transitionCounts, transitionMatrix] = ...
    empirical_transition_matrix(region, numberOfRegions);

predictedNextRegion = mod( ...
    region(1:end-1) - 1 + phaseSign, numberOfRegions) + 1;
observedNextRegion = region(2:end);
transitionAccuracy = mean( ...
    observedNextRegion(reliableStep) == ...
    predictedNextRegion(reliableStep));

fprintf('\nThree-region phase transport\n')
fprintf('----------------------------\n')
fprintf('Circular mean phase advance = %+10.6f rad\n', ...
    meanPhaseAdvance)
fprintf('Mismatch from selected target = %.6e rad\n', targetMismatch)
fprintf('Phase-residual concentration = %.6f\n', phaseConcentration)
fprintf('Reliable cyclic-transition fraction = %.6f\n', ...
    transitionAccuracy)
fprintf('Transition count matrix (rows: source; columns: destination):\n')
disp(transitionCounts)
fprintf('Row-normalized transition matrix:\n')
disp(transitionMatrix)

%% Spectral robustness across scalar observables

observableNames = {'x', 'dx/dt', 'x+dx/dt'};
observableData = [x(:, 1), x(:, 2), x(:, 1) + x(:, 2)];
observableSpectra = nan(numel(thetaGrid), numel(observableNames));
observablePeakLocation = nan(1, numel(observableNames));
observablePeakHeight = nan(1, numel(observableNames));

if runObservableRobustness
    fprintf('\nObservable robustness near 2*pi/3\n')
    fprintf('----------------------------------\n')

    for observableIndex = 1:numel(observableNames)
        standardizedObservable = standardize_observable( ...
            observableData(:, observableIndex));
        [PsiXObservable, PsiYObservable] = build_delay_dictionary( ...
            standardizedObservable, delayDimension);

        [~, density] = riggedDMD( ...
            PsiXObservable, PsiYObservable, uniformWeights, ...
            epsilon, thetaTarget, [], ...
            'order', kernelOrder, ...
            'g_coeffs', gCoefficients, ...
            'TH2', thetaGrid);

        density = real(squeeze(density));
        observableSpectra(:, observableIndex) = density(:);
        [observablePeakLocation(observableIndex), ...
         observablePeakHeight(observableIndex)] = local_maximum( ...
            thetaGrid, density, thetaTarget, 0.08*pi);

        fprintf('  %-8s: peak theta/pi = %.6f, height = %.6e\n', ...
            observableNames{observableIndex}, ...
            observablePeakLocation(observableIndex)/pi, ...
            observablePeakHeight(observableIndex))
    end
end

%% Spectral-density figure

spectralFigure = figure( ...
    'Name', 'Duffing Regime I spectral density', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);

plot(thetaGrid/pi, spectralDensity, ...
    'Color', [1 69/255 79/255], ...
    'LineWidth', 4)
hold on
xline(2/3, 'k--', '2\pi/3', 'LineWidth', 1.5)
xlabel('$\theta/\pi$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('Spectral density', ...
    'Interpreter', 'latex', 'FontSize', 20)
set(gca, 'FontSize', 16)
xlim([0 1])
grid on
box on

%% Observable-robustness figure

robustnessFigure = [];
if runObservableRobustness
    robustnessFigure = figure( ...
        'Name', 'Duffing Regime I observable robustness', ...
        'Color', 'w', ...
        'Position', [100 100 850 650]);
    hold on

    colours = lines(numel(observableNames));
    for observableIndex = 1:numel(observableNames)
        plot(thetaGrid/pi, ...
            observableSpectra(:, observableIndex), ...
            'LineWidth', 2, ...
            'Color', colours(observableIndex, :), ...
            'DisplayName', observableNames{observableIndex})
    end

    xline(2/3, 'k--', '2\pi/3', 'LineWidth', 1.5, ...
        'HandleVisibility', 'off')
    xlabel('$\theta/\pi$', ...
        'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
    ylabel('Spectral density', ...
        'Interpreter', 'latex', 'FontSize', 20)
    legend('Location', 'best', 'Interpreter', 'none')
    set(gca, 'FontSize', 16)
    xlim([0 1])
    grid on
    box on
end

%% Phase-partition figure

regionColours = [36 122 254; 230 159 0; 0 120 0]/255;

phaseFigure = figure( ...
    'Name', 'Duffing Regime I phase partition', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);
hold on

for regionIndex = 1:numberOfRegions
    inRegion = region == regionIndex;
    scatter(sectionData(inRegion, 1), sectionData(inRegion, 2), ...
        10, regionColours(regionIndex, :), 'filled', ...
        'DisplayName', sprintf('$R_%d$', regionIndex))
end

xlabel('$x$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\dot{x}$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'Interpreter', 'latex')
set(gca, 'FontSize', 16)
axis tight
grid on
box on

%% Wave-packet log-modulus figure

modulusFigure = figure( ...
    'Name', 'Duffing Regime I modulus', ...
    'Color', 'w', ...
    'Position', [100 100 850 650]);

scatter(sectionData(:, 1), sectionData(:, 2), ...
    12, logModulus, 'filled')
xlabel('$x$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\dot{x}$', ...
    'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
colourBar = colorbar;
colourBar.Label.String = '$\log|\phi_{\varepsilon,\theta}|$';
colourBar.Label.Interpreter = 'latex';
colourBar.Label.FontSize = 20;
colormap(parula)
set(gca, 'FontSize', 16)
axis tight
grid on
box on

%% Archive analysis results and metadata

analysis = struct();
analysis.generatedOn = datestr(now, 30);
analysis.observable = 'standardized position x';
analysis.observableMean = positionMean;
analysis.observableScale = positionScale;
analysis.N = N;
analysis.delayDimension = delayDimension;
analysis.nSnapshotPairs = nSnapshotPairs;
analysis.uniformWeight = uniformWeights(1);
analysis.uniformWeightSum = sum(uniformWeights);
analysis.epsilon = epsilon;
analysis.kernelOrder = kernelOrder;
analysis.thetaTarget = thetaTarget;
analysis.spectralPeakTheta = spectralPeakTheta;
analysis.spectralPeakValue = spectralPeakValue;
analysis.reliabilityQuantile = reliabilityQuantile;
analysis.modulusThreshold = modulusThreshold;
analysis.meanPhaseAdvance = meanPhaseAdvance;
analysis.phaseSign = phaseSign;
analysis.targetMismatch = targetMismatch;
analysis.phaseConcentration = phaseConcentration;
analysis.transitionCounts = transitionCounts;
analysis.transitionMatrix = transitionMatrix;
analysis.transitionAccuracy = transitionAccuracy;
analysis.riggedDMDImplementation = which('riggedDMD');

save(cacheFile, ...
    'x', 'sectionData', 'generationMetadata', 'analysis', ...
    'phi', 'phase', 'modulus', 'logModulus', ...
    'thetaGrid', 'spectralDensity', ...
    'region', 'delayDimension', 'epsilon', 'thetaTarget', ...
    'phaseSign', 'phaseResidual', 'reliableStep', ...
    'phaseConcentration', 'transitionCounts', 'transitionMatrix', ...
    'transitionAccuracy', 'observableNames', 'observableSpectra', ...
    'observablePeakLocation', 'observablePeakHeight', '-v7')

write_metadata_report(metadataFile, generationMetadata, analysis)

if saveFigures
    save_figure_pdf(spectralFigure, fullfile( ...
        outputDirectory, 'duffing_regime1_spectral_density.pdf'))
    if runObservableRobustness
        save_figure_pdf(robustnessFigure, fullfile( ...
            outputDirectory, ...
            'duffing_regime1_observable_robustness.pdf'))
    end
    save_figure_pdf(phaseFigure, fullfile( ...
        outputDirectory, 'duffing_regime1_phase_partition.pdf'))
    save_figure_pdf(modulusFigure, fullfile( ...
        outputDirectory, 'duffing_regime1_modulus.pdf'))
end

fprintf('\nSaved the section data and Regime I analysis to %s\n', cacheFile)
fprintf('Saved provenance and diagnostics to %s\n', metadataFile)
fprintf(['Run duffing_utils/duffing_upo_plot.m to regenerate the ' ...
    'six manuscript panels.\n'])

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
        'discardedCrossings', ...
        'integrationBlockLength', 'forcingPhase', 'forcingPeriod', ...
        'crossingRule', ...
        'eventLocationRule', 'retainedEventTimes', ...
        'maximumEventPhaseResidual', 'N'};

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
        metadata.forcingPeriod == config.forcingPeriod;

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
    positivePhaseDirection = config.omega > 0;

    assert(positivePhaseDirection, ...
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
    metadata.discardedTransientInterval = [ ...
        tStart, retainedEventTimes(1)];
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

function [standardized, center, scale] = standardize_observable(observable)
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

    fprintf(fileID, 'Forced Duffing Regime I stroboscopic archive\n');
    fprintf(fileID, 'Created: %s\n', metadata.createdOn);
    fprintf(fileID, 'Analyzed: %s\n\n', analysis.generatedOn);

    fprintf(fileID, ...
        'Parameters: alpha=%.16g, beta=%.16g, delta=%.16g, ', ...
        metadata.parameters.alpha, metadata.parameters.beta, ...
        metadata.parameters.delta);
    fprintf(fileID, 'omega=%.16g, gamma=%.16g\n', ...
        metadata.parameters.omega, metadata.parameters.gamma);
    fprintf(fileID, 'Initial condition [x,v,phi]: [%.16g, %.16g, %.16g]\n', ...
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
    fprintf(fileID, 'Maximum solver step: %.16g\n', metadata.maximumStep);
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
    fprintf(fileID, 'Smoothing epsilon: %.16g\n', analysis.epsilon);
    fprintf(fileID, 'Target theta: %.16g = %.16g*pi\n', ...
        analysis.thetaTarget, analysis.thetaTarget/pi);
    fprintf(fileID, 'Spectral peak theta: %.16g = %.16g*pi\n', ...
        analysis.spectralPeakTheta, analysis.spectralPeakTheta/pi);
    fprintf(fileID, 'Spectral peak value: %.16g\n', ...
        analysis.spectralPeakValue);
    fprintf(fileID, 'Circular mean phase advance: %.16g\n', ...
        analysis.meanPhaseAdvance);
    fprintf(fileID, 'Phase residual concentration: %.16g\n', ...
        analysis.phaseConcentration);
    fprintf(fileID, 'Reliable cyclic-transition fraction: %.16g\n\n', ...
        analysis.transitionAccuracy);

    fprintf(fileID, 'Transition counts (rows source, columns destination):\n');
    write_matrix(fileID, analysis.transitionCounts, '%12d');
    fprintf(fileID, 'Row-normalized transition matrix:\n');
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
