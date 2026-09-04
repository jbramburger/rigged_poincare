% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE KURAMOTO--SIVASHINSKY POINCARE MAP
% -------------------------------------------------------------------------
%
% This script reproduces the six panels in Figures KS_overview and
% KS_rigged for the K=32 Fourier--Galerkin approximation of the
% Kuramoto--Sivashinsky equation at nu=0.0298. It loads event-located
% Poincare-section data from ks_section_cache.mat; generate that archive
% once with generate_ks_section_data.m before running this analysis.
%
% The calculation proceeds as follows.
%
%   1. Load the archived upward crossings of a2=0. Each row of the section
%      matrix stores (a1,a3,...,a32) at a crossing with da2/dt>0. The cache
%      also records the crossing times, full crossing states, and numerical
%      provenance.
%
%   2. Retain sectionData(2:2:end,:), exactly as in the original analysis.
%      The first-return map P alternates between two symmetry-related
%      components, so these rows represent F=P^2 on one component. The
%      retained sample count is N=11544.
%
%   3. Compute the thin SVD of the uncentered N-by-31 section matrix.
%      The first two principal-component coordinates, X*v1 and X*v2, are
%      used to display the low-dimensional organization of the section.
%      No mean subtraction is applied before this SVD.
%
%   4. Construct d=20 delay-coordinate matrices from normalized PC1. The
%      uniform quadrature weights are normalized
%      over the N-d Hankel snapshot pairs; for N=11544, every weight is
%      1/11524.
%
%   5. Run riggedDMD twice on the same section data. The spectral density
%      and unwrapped phase use epsilon=0.5 and theta=pi/4. The log-modulus
%      uses a separate wave packet with epsilon=0.75 and theta=pi/4.
%
%   6. Order the section points with the original greedy nearest-neighbour
%      traversal, starting at the point with smallest PC1, and apply MATLAB
%      unwrap directly to the epsilon=0.5 packet phase in that order.
%      The unwrapped trace is retained in the phase plot; a robust median
%      trend and eight data-derived horizontal guide levels are overlaid to
%      make the plateau structure legible.
%
% The points A,...,F displayed in the section and log-modulus panels are
% loaded from per3asym.mat and projected into the same singular-coordinate
% basis as the chaotic section data. The point G is the preimage marker
% used in the paper. This script does not generate the separate periodic-
% orbit figures.
%
% A short integration is performed only for the attractor panel. The long
% event-data integration is performed separately and only once.
%
% REQUIREMENTS
%   - ks_section_cache.mat, created by generate_ks_section_data.m
%   - riggedDMD and its dependencies in main_routines/ or on the MATLAB path
%   - per3asym.mat in ks_upos/ or beside this script
%
% OUTPUT DATA
%   ks_section_metadata.txt        data-source and analysis diagnostics
%
% OUTPUT FIGURES
%   ks_a1a2.pdf               attractor and the section a2=0
%   ks_singular_values.pdf    normalized uncentered singular values
%   ks_psec.pdf               section data in the (PC1,PC2)-plane
%   ks_spec_measure.pdf       spectral density, epsilon=0.5
%   ks_unwrapped.pdf          packet phase under the curve ordering
%   ks_eigenfunction.pdf      wave-packet log-modulus, epsilon=0.75
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clearvars
close all
clc

%% Parameters

sectionCacheFileName = 'ks_section_cache.mat';
metadataFile = 'ks_section_metadata.txt';
outputDirectory = '.';

nu = 0.0298;
K = 32;

x0 = [ ...
    0;
    0.49744707;
    2.115777;
   -4.6417336;
   -3.915899;
   -5.449332;
    0.9944877;
    0.06847247;
   -0.20083837;
   -0.872909;
   -0.29952765;
   -0.13332021;
    0.06590602;
   -0.02766516;
    zeros(K-14,1)];

displayTimes = (0:1e-3:(50000-1)*1e-3).';
relativeTolerance = 1e-12;
absoluteTolerance = 1e-12*ones(1,K);

retainedParity = 2;
expectedArchivedRows = 23089;
expectedN = 11544;
expectedHankelRows = 11524;

delayDimension = 20;
spectralEpsilon = 0.6;
wavePacketEpsilon = 0.75;
wavePacketTheta = pi/4;
phaseGuideWindow = 301;
numberOfPhasePlateaus = 8;
thetaGrid = -pi:0.01:pi;
kernelOrder = 2;

if ~isfolder(outputDirectory)
    mkdir(outputDirectory)
end

scriptDirectory = fileparts(mfilename('fullpath'));
if isempty(scriptDirectory)
    scriptDirectory = pwd;
end

%% Locate the archived data and dependencies

add_riggeddmd_path(scriptDirectory)
sectionCacheFile = locate_file( ...
    sectionCacheFileName,scriptDirectory,'');
periodThreeFile = locate_file('per3asym.mat', scriptDirectory, 'ks_upos');

%% Load the Poincare-section data used for the published figures

loaded = load(sectionCacheFile,'cache');
assert(isfield(loaded,'cache'), ...
    'The section archive must contain a cache structure.')

sectionCache = loaded.cache;
requiredCacheFields = { ...
    'nu','modes','sectionData','crossingRule','firstReturnCount', ...
    'secondReturnCount','retainedCrossingTimes', ...
    'maximumSectionResidual','minimumCrossingSpeed'};
assert(all(isfield(sectionCache,requiredCacheFields)), ...
    'The section cache is missing required data or metadata.')
assert(sectionCache.modes == K && sectionCache.nu == nu, ...
    'The archived section data use different model parameters.')
assert(size(sectionCache.sectionData,1) == expectedArchivedRows && ...
    size(sectionCache.sectionData,2) == K-1, ...
    'The archived section matrix must have size 23089-by-31.')

archivedSectionData = sectionCache.sectionData;
X = archivedSectionData(retainedParity:2:end,:);
N = size(X,1);

assert(strcmp(sectionCache.crossingRule,'a2=0 with da2/dt>0'), ...
    'The supplied cache was not generated with the required crossing rule.')
assert(sectionCache.firstReturnCount == expectedArchivedRows, ...
    'The cache metadata and archived section matrix are inconsistent.')
assert(sectionCache.secondReturnCount == expectedN, ...
    'The cache metadata do not report N=11544 for F=P^2.')

assert(N == expectedN, ...
    'The retained F=P^2 data must contain N=11544 states.')

fprintf('Event-located section rows: %d\n',size(archivedSectionData,1))
fprintf('Retained section points for F=P^2: N = %d\n', N)
fprintf('Retained row rule: sectionData(2:2:end,:)\n')
fprintf('Crossing-time interval: [%.16g, %.16g]\n', ...
    sectionCache.retainedCrossingTimes(1), ...
    sectionCache.retainedCrossingTimes(end))
fprintf('Maximum retained |a2|: %.6e\n', ...
    sectionCache.maximumSectionResidual)
fprintf('Minimum retained da2/dt: %.16g\n', ...
    sectionCache.minimumCrossingSpeed)

%% Verify the alternating components and retained parity

expectedReducedData = archivedSectionData(retainedParity:2:end,:);
assert(isequal(X,expectedReducedData), ...
    'The reduced data do not match the documented retained parity.')

oddComponent = archivedSectionData(1:2:end,:);
evenComponent = archivedSectionData(2:2:end,:);
sectionIndices = [1 3:K];
symmetrySigns = (-1).^sectionIndices;
symmetryTransformedEven = evenComponent.*symmetrySigns;

oddCentroid = mean(oddComponent,1);
evenCentroid = mean(evenComponent,1);
transformedEvenCentroid = mean(symmetryTransformedEven,1);

componentMidpoint = 0.5*(oddCentroid+evenCentroid);
componentDirection = oddCentroid-evenCentroid;
componentScores = ...
    (archivedSectionData-componentMidpoint)*componentDirection.';
oddComponentPurity = mean(componentScores(1:2:end) > 0);
evenComponentPurity = mean(componentScores(2:2:end) < 0);
componentAlternationRate = mean( ...
    componentScores(1:end-1).*componentScores(2:end) < 0);

symmetryCentroidError = norm(oddCentroid-transformedEvenCentroid) ...
    / max(norm(oddCentroid),eps);
componentCentroidSeparation = norm(oddCentroid-evenCentroid);

assert(componentAlternationRate > 0.99 && ...
    oddComponentPurity > 0.99 && evenComponentPurity > 0.99, ...
    'The event-located first-return sequence does not alternate components.')

fprintf('Odd/even component sizes: %d and %d\n', ...
    size(oddComponent,1),size(evenComponent,1))
fprintf('Symmetry-transformed centroid error: %.6e\n', ...
    symmetryCentroidError)
fprintf('Odd/even centroid separation: %.6e\n', ...
    componentCentroidSeparation)
fprintf('Empirical component alternation rate: %.10f\n', ...
    componentAlternationRate)
fprintf('Odd/even component purities: %.10f and %.10f\n', ...
    oddComponentPurity,evenComponentPurity)
fprintf('Retained parity verified exactly.\n')

%% Plot the attractor and section in the (a1,a2)-plane

displayOptions = odeset( ...
    'RelTol',relativeTolerance, ...
    'AbsTol',absoluteTolerance);

[~,displayStates] = ode45( ...
    @(t,state) ks_rhs(state,nu,K), ...
    displayTimes,x0,displayOptions);

fig = figure('Color','w');
hold on

plot(displayStates(:,1),displayStates(:,2),'k','LineWidth',2)

a1Line = -1.5:0.01:1.5;
plot(a1Line, zeros(size(a1Line)), ...
    'Color', [0 120/255 0], 'LineWidth', 5)

xlabel('$a_1$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
ylabel('$a_2$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
set(gca,'FontSize',16)
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(fig, fullfile(outputDirectory,'ks_a1a2.pdf'))

%% Compute the principal-component coordinates

[U,S,V] = svd(X,'econ');
singularValues = diag(S);
principalCoordinates = U*S;

pc1 = principalCoordinates(:,1);
pc2 = principalCoordinates(:,2);

%% Plot the singular-value spectrum

fig = figure('Color','w');
singularWeights = singularValues/sum(singularValues);

bar(singularWeights, ...
    'FaceColor',[36/255 122/255 254/255], ...
    'EdgeColor','none')
hold on

smallValues = singularWeights < 1e-3;
plot(find(smallValues), singularWeights(smallValues), 'k*', ...
    'MarkerSize',8, 'LineWidth',1.5)

xlabel('Index', 'FontSize',16)
ylabel('Normalized singular values', 'FontSize',16)
set(gca,'FontSize',16)
grid on
box on

save_figure_pdf(fig, ...
    fullfile(outputDirectory,'ks_singular_values.pdf'))

%% Project the period-three section points

markers = load_period_three_markers(periodThreeFile,V);

fprintf('Projected section coordinates:\n')
fprintf('  A = [%.10f, %.10f]\n', markers.A)
fprintf('  B = [%.10f, %.10f]\n', markers.B)
fprintf('  C = [%.10f, %.10f]\n', markers.C)
fprintf('  D = [%.10f, %.10f]\n', markers.D)
fprintf('  E = [%.10f, %.10f]\n', markers.E)
fprintf('  F = [%.10f, %.10f]\n', markers.F)
fprintf('  G = [%.10f, %.10f]\n', markers.G)

%% Plot the reduced Poincare section

fig = figure('Color','w');
hold on

plot(pc1,pc2,'k.','MarkerSize',10)

plot(markers.DEF(:,1),markers.DEF(:,2),'diamond', ...
    'MarkerSize',16, ...
    'Color',[0 120/255 0], ...
    'MarkerFaceColor',[0 120/255 0])

plot(markers.CAB(:,1),markers.CAB(:,2),'square', ...
    'MarkerSize',16, ...
    'Color',[1 69/255 79/255], ...
    'MarkerFaceColor',[1 69/255 79/255])

plot(markers.G(1),markers.G(2),'pentagram', ...
    'MarkerSize',20, ...
    'Color',[151/255 1/255 200/255], ...
    'MarkerFaceColor',[151/255 14/255 200/255])

text(markers.C(1)-0.20,markers.C(2)+0.30,'C','FontSize',20)
text(markers.A(1)-0.20,markers.A(2)+0.30,'A','FontSize',20)
text(markers.B(1)-0.20,markers.B(2)+0.30,'B','FontSize',20)
text(markers.F(1)-0.15,markers.F(2)+0.27,'F','FontSize',20)
text(markers.D(1)-0.18,markers.D(2)+0.27,'D','FontSize',20)
text(markers.E(1)+0.05,markers.E(2)-0.27,'E','FontSize',20)
text(markers.G(1)+0.05,markers.G(2)-0.27,'G','FontSize',20)

xlabel('$\mathrm{PC}_1$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
ylabel('$\mathrm{PC}_2$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
set(gca,'FontSize',16)
grid on
box on

save_figure_pdf(fig, fullfile(outputDirectory,'ks_psec.pdf'))

%% Build the d=20 delay-coordinate matrices

observable = pc1-mean(pc1);
observable = observable/std(observable);

[PsiX,PsiY] = build_delay_dictionary(observable,delayDimension);
nSnapshotPairs = size(PsiX,1);
uniformWeights = ones(nSnapshotPairs,1)/nSnapshotPairs;

assert(nSnapshotPairs == N-delayDimension, ...
    'The number of Hankel rows must equal N-d.')
assert(nSnapshotPairs == expectedHankelRows, ...
    'The d=20 Hankel matrices must contain 11524 rows.')

fprintf('Hankel snapshot pairs N-d = %d\n', nSnapshotPairs)
fprintf('Uniform quadrature weight = 1/%d = %.16g\n', ...
    nSnapshotPairs,uniformWeights(1))
fprintf('Sum of weights = %.16g\n',sum(uniformWeights))

gCoefficients = zeros(delayDimension,1);
gCoefficients(1) = 1;

%% Spectral density: epsilon=0.5

fprintf('Running spectral density with epsilon = %.2f.\n', ...
    spectralEpsilon)

[gMODES,spectralDensity] = riggedDMD( ...
    PsiX,PsiY,uniformWeights,spectralEpsilon,wavePacketTheta,[], ...
    'order',kernelOrder, ...
    'g_coeffs',gCoefficients, ...
    'TH2',thetaGrid);

spectralDensity = real(squeeze(spectralDensity));
gMODES = squeeze(gMODES);
C = PsiX*gMODES;

positiveMask = thetaGrid >= 0 & thetaGrid <= pi;
positiveTheta = thetaGrid(positiveMask);
positiveDensity = spectralDensity(positiveMask);
[spectralPeakValue,peakIndex] = max(positiveDensity);
spectralPeakTheta = positiveTheta(peakIndex);

fprintf('Spectral peak: theta = %.10f = %.10f*pi\n', ...
    spectralPeakTheta,spectralPeakTheta/pi)

fig = figure('Color','w');
plot(thetaGrid,spectralDensity, ...
    'Color',[1 69/255 79/255], 'LineWidth',4)

xlabel('$\theta$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'FontSize',16)
xlim([0 pi])
grid on
box on

save_figure_pdf(fig, ...
    fullfile(outputDirectory,'ks_spec_measure.pdf'))

%% Wave packet: epsilon=0.75 and theta=pi/4

fprintf('Running wave packet with epsilon = %.2f and theta = pi/4.\n', ...
    wavePacketEpsilon)

[wavePacketCoefficients,~] = riggedDMD( ...
    PsiX,PsiY,uniformWeights,wavePacketEpsilon,wavePacketTheta,[], ...
    'order',kernelOrder, ...
    'g_coeffs',gCoefficients, ...
    'TH2',[]);

wavePacketCoefficients = squeeze(wavePacketCoefficients);
wavePacket = PsiX*wavePacketCoefficients;
wavePacketPhase = angle(wavePacket);
wavePacketLogModulus = log(abs(wavePacket));

%% Plot the spectral packet phase using the original ordering and unwrap construction

% This is the phase construction used in the supplied reference script:
% form the phase of C=PsiX*gMODES from the epsilon=0.5 run, wrap it only for
% the colour plot if needed, then order the section and unwrap the original
% angle values directly.
theta2 = angle(C);
theta2 = mod(theta2,2*pi);
P = [ ...
    pc1(1:length(theta2)), ...
    pc2(1:length(theta2))];
ord = greedy_curve_order(P);
theta2 = angle(C);
theta_ord = unwrap(theta2(ord));
phaseGuide = movmedian(theta_ord,phaseGuideWindow);
phasePlateauLevels = estimate_phase_plateau_levels( ...
    phaseGuide,numberOfPhasePlateaus);

fig = figure('Color','w');
plot(theta_ord,'Color',[170/255 205/255 1], 'LineWidth',0.5)
hold on
plot(phaseGuide,'Color',[36/255 122/255 254/255], 'LineWidth',1.5)
for levelIndex = 1:numel(phasePlateauLevels)
    yline(phasePlateauLevels(levelIndex),'--', ...
        'Color',[0.45 0.45 0.45], 'LineWidth',1)
end

xlabel('Ordered data point index', ...
    'FontSize',24, 'FontWeight','bold')
ylabel('Unwrapped phase', ...
    'FontSize',24, 'FontWeight','bold')
set(gca,'FontSize',16)
xlim([1 numel(theta_ord)])
grid on
box on

save_figure_pdf(fig, fullfile(outputDirectory,'ks_unwrapped.pdf'))

%% Plot the wave-packet log-modulus

pc1Plot = pc1(1:nSnapshotPairs);
[pc1Sorted,sortIndex] = sort(pc1Plot);
logModulusSorted = wavePacketLogModulus(sortIndex);

fig = figure('Color','w');
hold on

redColor    = [1 69/255 79/255];
greenColor  = [0 120/255 0];
purpleColor = [151/255 1/255 200/255];

xline(markers.A(1),'--','Color',redColor,'LineWidth',2)
xline(markers.B(1),'--','Color',redColor,'LineWidth',2)
xline(markers.C(1),'--','Color',redColor,'LineWidth',2)

xline(markers.D(1),'--','Color',greenColor,'LineWidth',2)
xline(markers.E(1),'--','Color',greenColor,'LineWidth',2)
xline(markers.F(1),'--','Color',greenColor,'LineWidth',2)

xline(markers.G(1),'--','Color',purpleColor,'LineWidth',2)

plot(pc1Sorted,logModulusSorted,'k-','LineWidth',1.5)

xlabel('$\mathrm{PC}_1$', 'Interpreter','latex', ...
    'FontSize',24, 'FontWeight','bold')
ylabel('$\log|\varphi_{\varepsilon,\theta}|$', ...
    'Interpreter','latex', 'FontSize',24, 'FontWeight','bold')
set(gca,'FontSize',16)
axis([7.75 10.6 -5.5 0])
grid on
box on

% Label the vertical lines just above the lower horizontal axis,
% slightly to the left of each line.
ax = gca;
xLimits = xlim(ax);
yLimits = ylim(ax);

labelDx = 0.010*diff(xLimits);
labelY  = yLimits(1) + 0.025*diff(yLimits);

text(markers.A(1)-labelDx,labelY,'A', ...
    'Color',redColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
text(markers.B(1)-labelDx,labelY,'B', ...
    'Color',redColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
text(markers.C(1)-labelDx,labelY,'C', ...
    'Color',redColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

text(markers.D(1)-labelDx,labelY,'D', ...
    'Color',greenColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
text(markers.E(1)-labelDx,labelY,'E', ...
    'Color',greenColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
text(markers.F(1)-labelDx,labelY,'F', ...
    'Color',greenColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

text(markers.G(1)-labelDx,labelY,'G', ...
    'Color',purpleColor,'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

save_figure_pdf(fig, ...
    fullfile(outputDirectory,'ks_eigenfunction.pdf'))

%% Write diagnostics

diagnostics = struct();
diagnostics.N = N;
diagnostics.symmetryCentroidError = symmetryCentroidError;
diagnostics.componentCentroidSeparation = componentCentroidSeparation;
diagnostics.componentAlternationRate = componentAlternationRate;
diagnostics.oddComponentPurity = oddComponentPurity;
diagnostics.evenComponentPurity = evenComponentPurity;
diagnostics.delayDimension = delayDimension;
diagnostics.nSnapshotPairs = nSnapshotPairs;
diagnostics.uniformWeight = uniformWeights(1);
diagnostics.spectralEpsilon = spectralEpsilon;
diagnostics.phaseEpsilon = spectralEpsilon;
diagnostics.spectralPeakTheta = spectralPeakTheta;
diagnostics.spectralPeakValue = spectralPeakValue;
diagnostics.wavePacketEpsilon = wavePacketEpsilon;
diagnostics.wavePacketTheta = wavePacketTheta;
diagnostics.phaseOrderStartIndex = ord(1);
diagnostics.phaseOrderStartPC1 = P(ord(1),1);
diagnostics.phaseOrderStartPC2 = P(ord(1),2);
diagnostics.phaseGuideWindow = phaseGuideWindow;
diagnostics.phasePlateauLevels = phasePlateauLevels;

dataInfo = struct();
dataInfo.sourceFile = sectionCacheFile;
dataInfo.nu = nu;
dataInfo.K = K;
dataInfo.archivedRows = size(archivedSectionData,1);
dataInfo.archivedColumns = size(archivedSectionData,2);
dataInfo.retainedParity = retainedParity;
dataInfo.cache = sectionCache;
dataInfo.displayInterval = [displayTimes(1),displayTimes(end)];
dataInfo.relativeTolerance = relativeTolerance;
dataInfo.absoluteTolerance = absoluteTolerance(1);

write_metadata(fullfile(outputDirectory,metadataFile),dataInfo,diagnostics)

fprintf('Finished writing the six KS figure panels.\n')

%% Local functions

function derivative = ks_rhs(state,nu,K)
    derivative = zeros(K,1);

    for k = 1:K
        derivative(k) = k^2*(1-nu*k^2)*state(k);

        for n = 1:(K-k)
            derivative(k) = derivative(k) ...
                + 0.5*k*state(n)*state(n+k);
        end

        for j = 1:(k-1)
            derivative(k) = derivative(k) ...
                - 0.25*k*state(j)*state(k-j);
        end
    end
end

function markers = load_period_three_markers(fileName,V)
    loaded = load(fileName);

    assert(isfield(loaded,'per3asymPsec1') && ...
        isfield(loaded,'per3asymPsec2'), ...
        ['The period-three file must contain per3asymPsec1 and ' ...
        'per3asymPsec2.'])

    orbit1 = loaded.per3asymPsec1;
    orbit2 = loaded.per3asymPsec2;

    analysisDimension = size(V,1);
    padded1 = zeros(size(orbit1,1),analysisDimension);
    padded2 = zeros(size(orbit2,1),analysisDimension);
    padded1(:,1:size(orbit1,2)) = orbit1;
    padded2(:,1:size(orbit2,2)) = orbit2;

    projected1 = padded1*V(:,1:2);
    projected2 = padded2*V(:,1:2);

    markers = struct();
    markers.DEF = projected1;
    markers.CAB = projected2;
    markers.F = projected1(1,:);
    markers.D = projected1(2,:);
    markers.E = projected1(3,:);
    markers.C = projected2(1,:);
    markers.A = projected2(2,:);
    markers.B = projected2(3,:);
    markers.G = [8.15 -1.98];
end

function [PsiX,PsiY] = build_delay_dictionary(observable,dimension)
    nRows = numel(observable)-dimension;
    rowOffsets = (1:nRows).';
    delayOffsets = 0:(dimension-1);
    indices = rowOffsets+delayOffsets;

    PsiX = observable(indices);
    PsiY = observable(indices+1);
end

function order = greedy_curve_order(points)
    numberOfPoints = size(points,1);
    order = zeros(numberOfPoints,1);
    used = false(numberOfPoints,1);

    % Match the original figure by beginning at the leftmost section point.
    [~,startIndex] = min(points(:,1));
    order(1) = startIndex;
    used(startIndex) = true;

    for j = 2:numberOfPoints
        lastPoint = points(order(j-1),:);
        distanceSquared = sum((points-lastPoint).^2,2);
        distanceSquared(used) = inf;
        [~,nextIndex] = min(distanceSquared);
        order(j) = nextIndex;
        used(nextIndex) = true;
    end
end

function levels = estimate_phase_plateau_levels(signal,numberOfLevels)
    % Estimate robust horizontal guide levels by one-dimensional k-means.
    signal = signal(:);
    signal = signal(isfinite(signal));
    assert(numel(signal) >= numberOfLevels, ...
        'The phase trace is too short to estimate plateau levels.')

    sortedSignal = sort(signal);
    initialIndices = round(linspace(1,numel(sortedSignal), ...
        numberOfLevels+2));
    levels = sortedSignal(initialIndices(2:end-1)).';

    for iteration = 1:100
        [~,labels] = min(abs(bsxfun(@minus,signal,levels)),[],2);
        updatedLevels = levels;
        for levelIndex = 1:numberOfLevels
            members = signal(labels == levelIndex);
            if ~isempty(members)
                updatedLevels(levelIndex) = median(members);
            end
        end

        if max(abs(updatedLevels-levels)) < 1e-10
            levels = updatedLevels;
            break
        end
        levels = updatedLevels;
    end

    levels = sort(levels);
end

function add_riggeddmd_path(scriptDirectory)
    candidates = { ...
        fullfile(pwd,'main_routines'), ...
        fullfile(scriptDirectory,'main_routines'), ...
        fullfile(scriptDirectory,'..','main_routines')};

    for j = 1:numel(candidates)
        if isfolder(candidates{j})
            addpath(candidates{j})
            break
        end
    end

    assert(exist('riggedDMD','file') == 2, ...
        'Could not find riggedDMD.m in main_routines or on the MATLAB path.')
end

function fileName = locate_file(baseName,scriptDirectory,subfolder)
    candidates = { ...
        fullfile(pwd,subfolder,baseName), ...
        fullfile(scriptDirectory,subfolder,baseName), ...
        fullfile(scriptDirectory,'..',subfolder,baseName), ...
        fullfile(pwd,baseName), ...
        fullfile(scriptDirectory,baseName)};

    fileName = '';
    for j = 1:numel(candidates)
        if isfile(candidates{j})
            fileName = candidates{j};
            break
        end
    end

    assert(~isempty(fileName), ['Could not locate ' baseName '.'])
end

function write_metadata(fileName,dataInfo,diagnostics)
    fileID = fopen(fileName,'w');
    assert(fileID >= 0,'Could not open the metadata file.')
    cleanup = onCleanup(@() fclose(fileID)); 

    fprintf(fileID,'Kuramoto--Sivashinsky Poincare-section analysis\n\n');
    fprintf(fileID,'Section-data source = %s\n',dataInfo.sourceFile);
    fprintf(fileID,'nu = %.16g\n',dataInfo.nu);
    fprintf(fileID,'K = %d\n',dataInfo.K);
    fprintf(fileID,'Archived matrix size = %d x %d\n', ...
        dataInfo.archivedRows,dataInfo.archivedColumns);
    fprintf(fileID,'Section-data solver = %s\n',dataInfo.cache.solver);
    fprintf(fileID,'Section relative tolerance = %.16g\n', ...
        dataInfo.cache.relativeTolerance);
    fprintf(fileID,'Section absolute tolerance = %.16g per component\n', ...
        dataInfo.cache.absoluteTolerance(1));
    fprintf(fileID,'Section maximum step = %.16g\n', ...
        dataInfo.cache.maximumStep);
    fprintf(fileID,'Requested section-data interval = [%.16g, %.16g]\n', ...
        dataInfo.cache.requestedIntegrationInterval);
    fprintf(fileID,'Actual section-data interval = [%.16g, %.16g]\n', ...
        dataInfo.cache.actualIntegrationInterval);
    fprintf(fileID,'Crossing rule = %s\n',dataInfo.cache.crossingRule);
    fprintf(fileID,'Event-location rule = %s\n', ...
        dataInfo.cache.eventLocationRule);
    fprintf(fileID,'Interpolation validation interval = [%.16g, %.16g]\n', ...
        dataInfo.cache.validationInterval);
    fprintf(fileID,'Validation crossings = %d\n', ...
        dataInfo.cache.validationCrossingCount);
    fprintf(fileID,'Maximum validation time error = %.16g\n', ...
        dataInfo.cache.maximumValidationTimeError);
    fprintf(fileID,'Maximum validation state error = %.16g\n', ...
        dataInfo.cache.maximumValidationStateError);
    fprintf(fileID,'Transient rule = %s\n',dataInfo.cache.transientRule);
    fprintf(fileID,'Discarded initial crossings = %d\n', ...
        dataInfo.cache.discardedInitialCrossings);
    fprintf(fileID, ...
        'Retained row rule for F=P^2 = sectionData(2:2:end,:)\n');
    fprintf(fileID,'Retained parity = %d\n',dataInfo.retainedParity);
    fprintf(fileID,'Retained count for F=P^2 = %d\n\n',diagnostics.N);

    fprintf(fileID,'Attractor-panel solver = ode45\n');
    fprintf(fileID,'Attractor-panel interval = [%.16g, %.16g]\n', ...
        dataInfo.displayInterval);
    fprintf(fileID,'Relative tolerance = %.16g\n', ...
        dataInfo.relativeTolerance);
    fprintf(fileID,'Absolute tolerance = %.16g in each component\n\n', ...
        dataInfo.absoluteTolerance);

    fprintf(fileID,'Symmetry-transformed centroid error = %.16g\n', ...
        diagnostics.symmetryCentroidError);
    fprintf(fileID,'Odd/even centroid separation = %.16g\n\n', ...
        diagnostics.componentCentroidSeparation);
    fprintf(fileID,'Component alternation rate = %.16g\n', ...
        diagnostics.componentAlternationRate);
    fprintf(fileID,'Odd component purity = %.16g\n', ...
        diagnostics.oddComponentPurity);
    fprintf(fileID,'Even component purity = %.16g\n\n', ...
        diagnostics.evenComponentPurity);

    fprintf(fileID,'Delay dimension = %d\n',diagnostics.delayDimension);
    fprintf(fileID,'Hankel rows = %d\n',diagnostics.nSnapshotPairs);
    fprintf(fileID,'Uniform weight = %.16g\n',diagnostics.uniformWeight);
    fprintf(fileID,'Density epsilon = %.16g\n', ...
        diagnostics.spectralEpsilon);
    fprintf(fileID,'Unwrapped-phase epsilon = %.16g\n', ...
        diagnostics.phaseEpsilon);
    fprintf(fileID,'Phase ordering = leftmost-start greedy nearest-neighbour\n');
    fprintf(fileID,'Phase-order start index = %d\n', ...
        diagnostics.phaseOrderStartIndex);
    fprintf(fileID,'Phase-order start (PC1,PC2) = [%.16g, %.16g]\n', ...
        diagnostics.phaseOrderStartPC1,diagnostics.phaseOrderStartPC2);
    fprintf(fileID,'Phase-guide median window = %d\n', ...
        diagnostics.phaseGuideWindow);
    fprintf(fileID,'Phase-guide plateau levels =');
    fprintf(fileID,' %.16g',diagnostics.phasePlateauLevels);
    fprintf(fileID,'\n');
    fprintf(fileID,'Spectral peak angle = %.16g\n', ...
        diagnostics.spectralPeakTheta);
    fprintf(fileID,'Spectral peak value = %.16g\n', ...
        diagnostics.spectralPeakValue);
    fprintf(fileID,'Log-modulus epsilon = %.16g\n', ...
        diagnostics.wavePacketEpsilon);
    fprintf(fileID,'Wave-packet theta = %.16g\n', ...
        diagnostics.wavePacketTheta);
end

function save_figure_pdf(figHandle,fileName)
    % Save a tightly sized vector PDF with no default paper margins.
    [fileDirectory,~,extension] = fileparts(fileName);
    if isempty(extension)
        fileName = [fileName '.pdf'];
    elseif ~strcmpi(extension,'.pdf')
        error('Figure output must have a .pdf extension.')
    end

    if ~isempty(fileDirectory) && ~isfolder(fileDirectory)
        mkdir(fileDirectory)
    end

    set(figHandle,'PaperUnits','centimeters')
    set(figHandle,'Units','centimeters')
    position = get(figHandle,'Position');
    set(figHandle,'PaperSize',[position(3) position(4)])
    set(figHandle,'PaperPositionMode','manual')
    set(figHandle,'PaperPosition',[0 0 position(3) position(4)])

    print(figHandle,fileName,'-dpdf','-painters')
end
