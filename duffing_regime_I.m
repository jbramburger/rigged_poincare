% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE FORCED DUFFING POINCARE MAP: REGIME I
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script analyzes the stroboscopic Poincare map of the periodically
% forced Duffing oscillator in Regime I.  It loads previously generated
% section data, constructs a delay-coordinate dictionary, and applies
% rigged Dynamic Mode Decomposition (riggedDMD) to identify the dominant
% three-fold organization of the return dynamics.
%
% The script computes:
%   1. The riggedDMD spectral density
%   2. A smoothed generalized eigenfunction near theta = 2*pi/3
%   3. The phase-based partition R1, R2, R3
%   4. The empirical transition matrix between phase regions
%   5. The phase-advance residual and its circular concentration
%   6. Spectral robustness across the observables x, xdot, and x+xdot
%
% MODEL
% The forced Duffing equation is
%
%   x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t),
%
% with
%   alpha = -1, beta = 1, delta = 0.3, omega = 1.2, gamma = 0.5.
%
% The Poincare map is obtained by sampling once per forcing period
% T = 2*pi/omega.  The supplied file
% duffing_utils/duffing_psec_I.mat must contain an N-by-2 array x whose
% columns are position and velocity.
%
% RIGGEDDMD IMPLEMENTATION
% The riggedDMD code used here is taken from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% REQUIRED FILES
%   duffing_utils/duffing_psec_I.mat
%                             : Poincare data stored as x = [position,velocity]
%   main_routines/            : riggedDMD and its dependencies
%
% OUTPUTS
%   duffing_regime1_results.mat
%   duffing_regime1_spectral_density.pdf
%   duffing_regime1_observable_robustness.pdf
%   duffing_regime1_phase_partition.pdf
%   duffing_regime1_modulus.pdf
%
% MAIN USER PARAMETERS
%   dim                       : number of delays
%   epsilon                   : riggedDMD smoothing parameter
%   order                     : rational-kernel order
%   thetaTarget               : generalized-eigenfunction frequency
%   reliabilityQuantile       : low-modulus fraction omitted from phase tests
%   runObservableRobustness   : compare three scalar observables
%   saveFigures               : save PDF figures
%
% NOTES
%   - The Koopman computations use only the sampled Poincare trajectory.
%   - The delay matrices are assembled directly, avoiding construction of
%     an unnecessary full Hankel matrix.
%   - The saved results file is used by duffing_utils/duffing_upo_plot.m.
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

dataFile = fullfile('duffing_utils','duffing_psec_I.mat');
routinesDirectory = 'main_routines';
outputDirectory = 'duffing_results';

dim = 20;
epsilon = 0.25;
order = 2;
thetaTarget = 2*pi/3;
thetaGrid = linspace(-pi,pi,1601);

reliabilityQuantile = 0.05;
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
assert(isfield(S,'x'),'duffing_psec_I.mat must contain the variable x.')
x = S.x;
assert(isnumeric(x) && size(x,2) == 2, ...
    'x must be an N-by-2 array [position, velocity].')
assert(all(isfinite(x(:))),'The Poincare data contain NaN or Inf values.')

%% Build the delay-coordinate dictionary

% The velocity observable is used for the phase partition shown in the
% paper.  The position and combined observables are checked below.
observable = normalize_observable(x(:,2));
[PsiX,PsiY] = delay_dictionary(observable,dim);

nSnapshots = size(PsiX,1);
sectionData = x(1:nSnapshots,:);
quadratureWeight = 1/size(x,1);

gCoeffs = zeros(dim,1);
gCoeffs(1) = 1;

%% Apply riggedDMD

[gMode,spectralDensity] = run_rigged_dmd( ...
    PsiX,PsiY,quadratureWeight,epsilon,thetaTarget,order,gCoeffs,thetaGrid);

phi = PsiX*gMode;
phase = angle(phi);
modulus = abs(phi);
logModulus = log(max(modulus,realmin));

%% Test phase advance and construct the three-region partition

phaseIncrement = wrap_to_pi(phase(2:end)-phase(1:end-1));
modulusThreshold = quantile(modulus,reliabilityQuantile);
reliableStep = modulus(1:end-1) > modulusThreshold & ...
               modulus(2:end) > modulusThreshold;

meanPhaseAdvance = angle(mean(exp(1i*phaseIncrement(reliableStep))));
targetMismatchPlus = abs(wrap_to_pi(meanPhaseAdvance-thetaTarget));
targetMismatchMinus = abs(wrap_to_pi(meanPhaseAdvance+thetaTarget));

if targetMismatchPlus <= targetMismatchMinus
    phaseSign = 1;
    phaseResidual = wrap_to_pi(phaseIncrement-thetaTarget);
    targetMismatch = targetMismatchPlus;
else
    phaseSign = -1;
    phaseResidual = wrap_to_pi(phaseIncrement+thetaTarget);
    targetMismatch = targetMismatchMinus;
end

reliableResidual = phaseResidual(reliableStep);
phaseConcentration = abs(mean(exp(1i*reliableResidual)));

numberOfRegions = 3;
regionWidth = 2*pi/numberOfRegions;
region = floor((phase+pi)/regionWidth)+1;
region = min(max(region,1),numberOfRegions);

regionShift = phaseSign;
predictedNext = mod(region(1:end-1)-1+regionShift,numberOfRegions)+1;
observedNext = region(2:end);

[transitionCounts,transitionMatrix] = ...
    empirical_transition_matrix(region,numberOfRegions);
transitionAccuracy = mean(observedNext(reliableStep)==predictedNext(reliableStep));

fprintf('\nDuffing Regime I riggedDMD analysis\n')
fprintf('-----------------------------------\n')
fprintf('Poincare points:                     %d\n',size(x,1))
fprintf('Delay dimension:                     %d\n',dim)
fprintf('Smoothing parameter:                 %.4f\n',epsilon)
fprintf('Circular mean phase advance:         %+10.6f rad\n',meanPhaseAdvance)
fprintf('Mismatch from selected target:       %.6e rad\n',targetMismatch)
fprintf('Phase-residual concentration:        %.6f\n',phaseConcentration)
fprintf('Reliable transition accuracy:        %.4f\n',transitionAccuracy)
disp('Transition count matrix:')
disp(transitionCounts)
disp('Row-normalized transition matrix:')
disp(transitionMatrix)

%% Check spectral robustness across observables

observableNames = {'x','dx/dt','x+dx/dt'};
observableData = [x(:,1),x(:,2),x(:,1)+x(:,2)];
observableSpectra = nan(numel(thetaGrid),numel(observableNames));
observablePeakLocation = nan(1,numel(observableNames));
observablePeakHeight = nan(1,numel(observableNames));

if runObservableRobustness
    for j = 1:numel(observableNames)
        gj = normalize_observable(observableData(:,j));
        [PsiXj,PsiYj] = delay_dictionary(gj,dim);
        [~,xij] = run_rigged_dmd(PsiXj,PsiYj,quadratureWeight,epsilon, ...
            thetaTarget,order,gCoeffs,thetaGrid);
        observableSpectra(:,j) = real(xij(:));
        [observablePeakLocation(j),observablePeakHeight(j)] = ...
            local_maximum(thetaGrid,observableSpectra(:,j),thetaTarget,0.08*pi);
    end

    fprintf('\nObservable robustness near 2*pi/3\n')
    for j = 1:numel(observableNames)
        fprintf('  %-8s: peak theta/pi = %.5f, height = %.4e\n', ...
            observableNames{j},observablePeakLocation(j)/pi,observablePeakHeight(j))
    end
end

%% Plot spectral density

figure('Name','Duffing Regime I spectral density','Color','w')
plot(thetaGrid/pi,real(spectralDensity),'Color',[1 69/255 79/255], ...
    'LineWidth',4)
xline(2/3,'k--','2\pi/3','LineWidth',1.5)
xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('Spectral density','Interpreter','latex','FontSize',20)
set(gca,'FontSize',16)
xlim([0 1])
grid on
box on

%% Plot observable robustness

if runObservableRobustness
    figure('Name','Duffing Regime I observable robustness','Color','w')
    hold on
    colours = lines(numel(observableNames));
    for j = 1:numel(observableNames)
        plot(thetaGrid/pi,observableSpectra(:,j),'LineWidth',2, ...
            'Color',colours(j,:),'DisplayName',observableNames{j})
    end
    xline(2/3,'k--','2\pi/3','LineWidth',1.5,'HandleVisibility','off')
    xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    ylabel('Spectral density','Interpreter','latex','FontSize',20)
    set(gca,'FontSize',16)
    xlim([0 1])
    legend('Location','best')
    grid on
    box on
end

%% Plot phase partition and generalized-eigenfunction modulus

regionColours = [36 122 254; 230 159 0; 0 120 0]/255;

figure('Name','Duffing Regime I phase partition','Color','w')
hold on
for j = 1:numberOfRegions
    idx = region == j;
    scatter(sectionData(idx,1),sectionData(idx,2),10,regionColours(j,:),'filled', ...
        'DisplayName',sprintf('$R_%d$',j))
end
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\dot{x}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
legend('Location','best','Interpreter','latex')
axis tight
grid on
box on

figure('Name','Duffing Regime I modulus','Color','w')
scatter(sectionData(:,1),sectionData(:,2),12,logModulus,'filled')
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\dot{x}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
colormap(parula)
colorbar
axis tight
grid on
box on

%% Save results and figures

outputPath = fullfile(scriptDirectory,outputDirectory);
if ~isfolder(outputPath)
    mkdir(outputPath)
end

resultsFile = fullfile(outputPath,'duffing_regime1_results.mat');
save(resultsFile,'x','sectionData','phi','phase','modulus','logModulus', ...
    'thetaGrid','spectralDensity','meanPhaseAdvance','targetMismatch', ...
    'region','dim','epsilon','thetaTarget','phaseSign','phaseResidual', ...
    'reliableStep','phaseConcentration','transitionCounts', ...
    'transitionMatrix','transitionAccuracy','observableNames', ...
    'observablePeakLocation','observablePeakHeight')

if saveFigures
    save_named_figure('Duffing Regime I spectral density',outputPath, ...
        'duffing_regime1_spectral_density')
    if runObservableRobustness
        save_named_figure('Duffing Regime I observable robustness',outputPath, ...
            'duffing_regime1_observable_robustness')
    end
    save_named_figure('Duffing Regime I phase partition',outputPath, ...
        'duffing_regime1_phase_partition')
    save_named_figure('Duffing Regime I modulus',outputPath, ...
        'duffing_regime1_modulus')
end

fprintf('\nSaved Regime I results to %s\n',resultsFile)

%% Local functions

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
