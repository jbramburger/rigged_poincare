% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE FORCED DUFFING POINCARE MAP: REGIME II
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script analyzes the stroboscopic Poincare map of the periodically
% forced Duffing oscillator in Regime II.  The dominant riggedDMD spectral
% component lies near theta = 6*pi/7 and induces a nonlocal seven-region
% transport cycle.  The script reproduces the smoothing-robustness table
% reported in the accompanying paper.
%
% The script computes:
%   1. Spectral densities over a prescribed smoothing-parameter sweep
%   2. The dominant peak near 6*pi/7 and weaker companion bumps
%   3. Smoothed generalized eigenfunctions at theta = 6*pi/7
%   4. The seven-region phase partition
%   5. The transport rule R_j -> R_{j+3} (mod 7)
%   6. Phase concentration, transition accuracy, and seven-step return rate
%   7. Enrichment of incorrect transitions near the low-modulus skeleton
%   8. Spectral robustness across x, xdot, and x+xdot observables
%
% MODEL
% The forced Duffing equation is
%
%   x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t),
%
% with
%   alpha = -1, beta = 0.25, delta = 0.1, omega = 2, gamma = 2.5.
%
% The Poincare map is obtained by sampling once per forcing period
% T = 2*pi/omega.  The supplied file
% duffing_utils/duffing_psec_II.mat must contain an N-by-2 array x whose
% columns are position and velocity.
%
% RIGGEDDMD IMPLEMENTATION
% The riggedDMD code used here is taken from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% REQUIRED FILES
%   duffing_utils/duffing_psec_II.mat
%                             : Poincare data stored as x = [position,velocity]
%   main_routines/            : riggedDMD and its dependencies
%
% OUTPUTS
%   duffing_regime2_results.mat
%   duffing_regime2_smoothing_sweep.pdf
%   duffing_regime2_observable_robustness.pdf
%   duffing_regime2_phase_partition.pdf
%   duffing_regime2_modulus.pdf
%   duffing_regime2_transition_matrix.pdf
%
% MAIN USER PARAMETERS
%   dim                       : number of delays
%   epsilonSweep              : smoothing values used in robustness tests
%   selectedEpsilon           : smoothing used in displayed phase plots
%   order                     : rational-kernel order
%   thetaTarget               : dominant frequency, here 6*pi/7
%   reliabilityQuantile       : lowest-modulus fraction omitted from phase tests
%   lowModulusQuantile        : retained low-modulus fraction used in enrichment
%   runObservableRobustness   : compare three scalar observables
%   saveFigures               : save PDF figures
%
% NOTES
%   - The reported transition-error enrichment is dynamically defined; it
%     does not use an arbitrary phase-residual percentile.
%   - A transition is reliable only when both endpoint moduli exceed the
%     prescribed reliability quantile.  The low-modulus decile is then
%     computed within this retained set.
%   - The delay matrices are assembled directly rather than through a full
%     N-by-N Hankel matrix.
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

dataFile = fullfile('duffing_utils','duffing_psec_II.mat');
routinesDirectory = 'main_routines';
outputDirectory = 'duffing_results';

dim = 20;
epsilonSweep = [0.05 0.10 0.15 0.20 0.25];
selectedEpsilon = 0.10;
order = 2;
thetaTarget = 6*pi/7;
thetaGrid = linspace(-pi,pi,1601);

reliabilityQuantile = 0.05;
lowModulusQuantile = 0.10;
runObservableRobustness = true;
saveFigures = true;

%% Locate dependencies and load Poincare data

scriptDirectory = fileparts(mfilename('fullpath'));
if isempty(scriptDirectory)
    scriptDirectory = pwd;
end

dataPath = resolve_path(scriptDirectory,dataFile,'file');
routinesPath = resolve_path(scriptDirectory,routinesDirectory,'directory');
addpath(routinesPath)
assert(exist('riggedDMD','file') == 2, ...
    'riggedDMD.m was not found in %s.',routinesPath)

S = load(dataPath);
assert(isfield(S,'x'),'duffing_psec_II.mat must contain the variable x.')
x = S.x;
assert(isnumeric(x) && size(x,2) == 2, ...
    'x must be an N-by-2 array [position, velocity].')
assert(all(isfinite(x(:))),'The Poincare data contain NaN or Inf values.')

%% Build delay-coordinate dictionary from the position observable

observable = normalize_observable(x(:,1));
[PsiX,PsiY] = delay_dictionary(observable,dim);

nSnapshots = size(PsiX,1);
sectionData = x(1:nSnapshots,:);
quadratureWeight = 1/size(x,1);

gCoeffs = zeros(dim,1);
gCoeffs(1) = 1;

%% Sweep the smoothing parameter

nEpsilon = numel(epsilonSweep);
spectralDensitySweep = nan(numel(thetaGrid),nEpsilon);
peakLocation = nan(nEpsilon,1);
peakHeight = nan(nEpsilon,1);
meanPhaseAdvance = nan(nEpsilon,1);
targetMismatch = nan(nEpsilon,1);
phaseConcentration = nan(nEpsilon,1);
medianAbsoluteResidual = nan(nEpsilon,1);
transitionAccuracy = nan(nEpsilon,1);
sevenStepReturn = nan(nEpsilon,1);
baselineTransportError = nan(nEpsilon,1);
lowModulusTransportError = nan(nEpsilon,1);
transportErrorEnrichment = nan(nEpsilon,1);
metricData = cell(nEpsilon,1);

for j = 1:nEpsilon
    epsilon = epsilonSweep(j);
    [gMode,xi] = run_rigged_dmd(PsiX,PsiY,quadratureWeight,epsilon, ...
        thetaTarget,order,gCoeffs,thetaGrid);

    phi = PsiX*gMode;
    metrics = transport_metrics(phi,thetaTarget,reliabilityQuantile, ...
        lowModulusQuantile,7);

    spectralDensitySweep(:,j) = real(xi(:));
    [peakLocation(j),peakHeight(j)] = local_maximum( ...
        thetaGrid,spectralDensitySweep(:,j),thetaTarget,0.06*pi);
    meanPhaseAdvance(j) = metrics.meanPhaseAdvance;
    targetMismatch(j) = metrics.targetMismatch;
    phaseConcentration(j) = metrics.phaseConcentration;
    medianAbsoluteResidual(j) = metrics.medianAbsoluteResidual;
    transitionAccuracy(j) = metrics.transitionAccuracy;
    sevenStepReturn(j) = metrics.sevenStepReturn;
    baselineTransportError(j) = metrics.baselineTransportError;
    lowModulusTransportError(j) = metrics.lowModulusTransportError;
    transportErrorEnrichment(j) = metrics.transportErrorEnrichment;
    metricData{j} = metrics;
end

selectedIndex = find(abs(epsilonSweep-selectedEpsilon) < 10*eps,1);
assert(~isempty(selectedIndex),'selectedEpsilon must be a member of epsilonSweep.')
selectedMetrics = metricData{selectedIndex};

robustnessTable = table(epsilonSweep(:),peakLocation/pi,phaseConcentration, ...
    medianAbsoluteResidual,transitionAccuracy,sevenStepReturn, ...
    transportErrorEnrichment,baselineTransportError,lowModulusTransportError, ...
    'VariableNames',{'epsilon','thetaMaxOverPi','phaseConcentration', ...
    'medianAbsResidual','pPlus3','p7','errorEnrichment', ...
    'baselineError','lowModulusError'});

fprintf('\nDuffing Regime II sevenfold transport diagnostics\n')
fprintf('-------------------------------------------------\n')
fprintf('Poincare points:                %d\n',size(x,1))
fprintf('Delay dimension:                %d\n',dim)
disp(robustnessTable)

%% Report companion spectral bumps

companionTargets = [2*pi/7 4*pi/7 6*pi/7];
companionLocation = nan(nEpsilon,3);
companionHeight = nan(nEpsilon,3);

fprintf('Local spectral maxima near 2*pi/7, 4*pi/7, and 6*pi/7:\n')
for j = 1:nEpsilon
    for k = 1:3
        [companionLocation(j,k),companionHeight(j,k)] = local_maximum( ...
            thetaGrid,spectralDensitySweep(:,j),companionTargets(k),0.06*pi);
    end
    fprintf('  epsilon = %.2f: theta/pi = [%7.4f %7.4f %7.4f]\n', ...
        epsilonSweep(j),companionLocation(j,:)/pi)
end

%% Check spectral robustness across observables at selected epsilon

observableNames = {'x','dx/dt','x+dx/dt'};
observableData = [x(:,1),x(:,2),x(:,1)+x(:,2)];
observableSpectra = nan(numel(thetaGrid),numel(observableNames));
observablePeakLocation = nan(1,numel(observableNames));
observablePeakHeight = nan(1,numel(observableNames));

if runObservableRobustness
    for j = 1:numel(observableNames)
        gj = normalize_observable(observableData(:,j));
        [PsiXj,PsiYj] = delay_dictionary(gj,dim);
        [~,xij] = run_rigged_dmd(PsiXj,PsiYj,quadratureWeight, ...
            selectedEpsilon,thetaTarget,order,gCoeffs,thetaGrid);
        observableSpectra(:,j) = real(xij(:));
        [observablePeakLocation(j),observablePeakHeight(j)] = local_maximum( ...
            thetaGrid,observableSpectra(:,j),thetaTarget,0.06*pi);
    end

    fprintf('\nObservable robustness near 6*pi/7 at epsilon = %.2f\n',selectedEpsilon)
    for j = 1:numel(observableNames)
        fprintf('  %-8s: peak theta/pi = %.5f, height = %.4e\n', ...
            observableNames{j},observablePeakLocation(j)/pi,observablePeakHeight(j))
    end
end

%% Plot smoothing sweep

figure('Name','Duffing Regime II smoothing sweep','Color','w')
hold on
colours = lines(nEpsilon);
for j = 1:nEpsilon
    plot(thetaGrid/pi,spectralDensitySweep(:,j),'LineWidth',1.8, ...
        'Color',colours(j,:),'DisplayName',sprintf('$\\varepsilon=%.2f$',epsilonSweep(j)))
end
xline(2/7,'k:','2\pi/7','HandleVisibility','off')
xline(4/7,'k:','4\pi/7','HandleVisibility','off')
xline(6/7,'k--','6\pi/7','LineWidth',1.5,'HandleVisibility','off')
xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('Spectral density','Interpreter','latex','FontSize',20)
set(gca,'FontSize',16)
xlim([0 1])
legend('Location','best','Interpreter','latex')
grid on
box on

%% Plot observable robustness

if runObservableRobustness
    figure('Name','Duffing Regime II observable robustness','Color','w')
    hold on
    observableColours = lines(numel(observableNames));
    for j = 1:numel(observableNames)
        plot(thetaGrid/pi,observableSpectra(:,j),'LineWidth',2, ...
            'Color',observableColours(j,:),'DisplayName',observableNames{j})
    end
    xline(6/7,'k--','6\pi/7','LineWidth',1.5,'HandleVisibility','off')
    xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    ylabel('Spectral density','Interpreter','latex','FontSize',20)
    set(gca,'FontSize',16)
    xlim([0 1])
    legend('Location','best')
    grid on
    box on
end

%% Plot selected phase partition, modulus, and transition matrix

phase = selectedMetrics.phase;
modulus = selectedMetrics.modulus;
region = selectedMetrics.region;
transitionMatrix = selectedMetrics.transitionMatrix;

figure('Name','Duffing Regime II phase partition','Color','w')
scatter(sectionData(:,1),sectionData(:,2),10,region,'filled')
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\dot{x}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
colormap(lines(7))
cb = colorbar;
cb.Ticks = 1:7;
axis tight
grid on
box on

figure('Name','Duffing Regime II modulus','Color','w')
scatter(sectionData(:,1),sectionData(:,2),12,log(max(modulus,realmin)),'filled')
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\dot{x}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
colormap(parula)
colorbar
axis tight
grid on
box on

figure('Name','Duffing Regime II transition matrix','Color','w')
imagesc(transitionMatrix,[0 1])
axis square
xlabel('Destination region','FontSize',16)
ylabel('Source region','FontSize',16)
title(sprintf('$R_j\\to R_{j%+d}\\; (\\mathrm{mod}\\;7)$', ...
    selectedMetrics.regionShift),'Interpreter','latex')
set(gca,'FontSize',16,'XTick',1:7,'YTick',1:7)
colorbar
hold on
for j = 1:7
    destination = mod(j-1+selectedMetrics.regionShift,7)+1;
    plot(destination,j,'ks','MarkerSize',14,'LineWidth',1.5)
end

%% Save results and figures

outputPath = fullfile(scriptDirectory,outputDirectory);
if ~isfolder(outputPath)
    mkdir(outputPath)
end

resultsFile = fullfile(outputPath,'duffing_regime2_results.mat');
save(resultsFile,'x','sectionData','dim','thetaTarget','thetaGrid', ...
    'epsilonSweep','selectedEpsilon','spectralDensitySweep','peakLocation', ...
    'peakHeight','companionLocation','companionHeight','robustnessTable', ...
    'metricData','selectedMetrics','observableNames','observablePeakLocation', ...
    'observablePeakHeight')

if saveFigures
    save_named_figure('Duffing Regime II smoothing sweep',outputPath, ...
        'duffing_regime2_smoothing_sweep')
    if runObservableRobustness
        save_named_figure('Duffing Regime II observable robustness',outputPath, ...
            'duffing_regime2_observable_robustness')
    end
    save_named_figure('Duffing Regime II phase partition',outputPath, ...
        'duffing_regime2_phase_partition')
    save_named_figure('Duffing Regime II modulus',outputPath, ...
        'duffing_regime2_modulus')
    save_named_figure('Duffing Regime II transition matrix',outputPath, ...
        'duffing_regime2_transition_matrix')
end

fprintf('\nSaved Regime II results to %s\n',resultsFile)

%% Local functions

function metrics = transport_metrics(phi,thetaTarget,reliabilityQuantile,lowModulusQuantile,nRegions)
    phase = angle(phi);
    modulus = abs(phi);
    phaseIncrement = wrap_to_pi(phase(2:end)-phase(1:end-1));

    modulusThreshold = quantile(modulus,reliabilityQuantile);
    reliableStep = modulus(1:end-1) > modulusThreshold & ...
                   modulus(2:end) > modulusThreshold;

    meanAdvance = angle(mean(exp(1i*phaseIncrement(reliableStep))));
    mismatchPlus = abs(wrap_to_pi(meanAdvance-thetaTarget));
    mismatchMinus = abs(wrap_to_pi(meanAdvance+thetaTarget));

    if mismatchPlus <= mismatchMinus
        phaseSign = 1;
        residual = wrap_to_pi(phaseIncrement-thetaTarget);
        mismatch = mismatchPlus;
    else
        phaseSign = -1;
        residual = wrap_to_pi(phaseIncrement+thetaTarget);
        mismatch = mismatchMinus;
    end

    reliableResidual = residual(reliableStep);
    phaseConcentration = abs(mean(exp(1i*reliableResidual)));

    regionWidth = 2*pi/nRegions;
    region = floor((phase+pi)/regionWidth)+1;
    region = min(max(region,1),nRegions);
    regionShift = round(phaseSign*thetaTarget/regionWidth);

    predictedNext = mod(region(1:end-1)-1+regionShift,nRegions)+1;
    observedNext = region(2:end);
    [counts,P] = empirical_transition_matrix(region,nRegions);

    transitionError = reliableStep & (observedNext ~= predictedNext);
    transitionAccuracy = mean(observedNext(reliableStep)==predictedNext(reliableStep));
    sevenStepReturn = mean(region(1:end-nRegions)==region(1+nRegions:end));

    stepModulus = min(modulus(1:end-1),modulus(2:end));
    lowThreshold = quantile(stepModulus(reliableStep),lowModulusQuantile);
    lowModulusStep = reliableStep & stepModulus <= lowThreshold;

    baselineError = sum(transitionError)/sum(reliableStep);
    conditionalError = sum(transitionError & lowModulusStep)/sum(lowModulusStep);
    errorEnrichment = conditionalError/baselineError;

    metrics = struct;
    metrics.phi = phi;
    metrics.phase = phase;
    metrics.modulus = modulus;
    metrics.phaseSign = phaseSign;
    metrics.phaseResidual = residual;
    metrics.reliableStep = reliableStep;
    metrics.meanPhaseAdvance = meanAdvance;
    metrics.targetMismatch = mismatch;
    metrics.phaseConcentration = phaseConcentration;
    metrics.medianAbsoluteResidual = median(abs(reliableResidual));
    metrics.region = region;
    metrics.regionShift = regionShift;
    metrics.transitionCounts = counts;
    metrics.transitionMatrix = P;
    metrics.transitionAccuracy = transitionAccuracy;
    metrics.sevenStepReturn = sevenStepReturn;
    metrics.transportError = transitionError;
    metrics.baselineTransportError = baselineError;
    metrics.lowModulusTransportError = conditionalError;
    metrics.transportErrorEnrichment = errorEnrichment;
end

function g = normalize_observable(g)
    g = g(:);
    g = (g-mean(g))/std(g);
end

function [PsiX,PsiY] = delay_dictionary(g,dim)
    g = g(:);
    nRows = numel(g)-dim;
    assert(nRows > 0,'The observable must contain more than dim samples.')
    PsiX = zeros(nRows,dim);
    PsiY = zeros(nRows,dim);
    for j = 1:dim
        PsiX(:,j) = g(j:j+nRows-1);
        PsiY(:,j) = g(j+1:j+nRows);
    end
end

function [mode,xi] = run_rigged_dmd(PsiX,PsiY,weight,epsilon,theta,order,gCoeffs,thetaGrid)
    [mode,xi] = riggedDMD(PsiX,PsiY,weight,epsilon,theta,[], ...
        'order',order,'g_coeffs',gCoeffs,'TH2',thetaGrid);
    mode = squeeze(mode);
end

function y = wrap_to_pi(x)
    y = mod(x+pi,2*pi)-pi;
end

function [counts,P] = empirical_transition_matrix(region,nRegions)
    counts = zeros(nRegions);
    for n = 1:numel(region)-1
        counts(region(n),region(n+1)) = counts(region(n),region(n+1))+1;
    end
    rowTotals = sum(counts,2);
    P = zeros(size(counts));
    occupied = rowTotals > 0;
    P(occupied,:) = counts(occupied,:)./rowTotals(occupied);
end

function [location,height] = local_maximum(thetaGrid,density,target,halfWidth)
    idx = thetaGrid >= target-halfWidth & thetaGrid <= target+halfWidth;
    thetaLocal = thetaGrid(idx);
    densityLocal = density(idx);
    [height,j] = max(densityLocal);
    location = thetaLocal(j);
end

function pathOut = resolve_path(scriptDirectory,name,pathType)
    pathOut = fullfile(scriptDirectory,name);
    if strcmp(pathType,'file') && ~isfile(pathOut)
        pathOut = name;
    elseif strcmp(pathType,'directory') && ~isfolder(pathOut)
        pathOut = name;
    end
    if strcmp(pathType,'file')
        assert(isfile(pathOut),'Could not find required file %s.',name)
    else
        assert(isfolder(pathOut),'Could not find required directory %s.',name)
    end
end

function save_named_figure(figureName,outputPath,fileName)
    fig = findobj(groot,'Type','figure','Name',figureName);
    assert(~isempty(fig),'Could not find figure named %s.',figureName)
    save_figure_pdf(fig(1),fullfile(outputPath,fileName))
end

function save_figure_pdf(figHandle,fileName)
    set(figHandle,'PaperUnits','centimeters')
    set(figHandle,'Units','centimeters')
    pos = get(figHandle,'Position');
    set(figHandle,'PaperSize',[pos(3) pos(4)])
    set(figHandle,'PaperPositionMode','manual')
    set(figHandle,'PaperPosition',[0 0 pos(3) pos(4)])
    print(figHandle,'-dpdf',fileName)
end
