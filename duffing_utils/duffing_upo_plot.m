% -------------------------------------------------------------------------
%  DUFFING_UPO_PLOT  Overlay Regime I unstable periodic orbits on the
%  riggedDMD phase and modulus fields.
% -------------------------------------------------------------------------
%
%  DESCRIPTION
%  This script combines the generalized eigenfunction computed by
%  duffing_regime_I.m with the periodic-orbit data computed by
%  duffing_upo_search.m.  For each requested minimal period it produces a
%  phase-partition plot and a log-modulus plot, with every stroboscopic
%  intersection of the corresponding unstable periodic orbits overlaid.
%
%  WORKFLOW
%    1. Load duffing_section_I_cache.mat from the repository root.
%    2. Load the previously computed duffing_upo_data.mat archive.
%    3. Select all orbit cycles having each requested minimal period.
%    4. Overlay their section points on the three-region phase partition.
%    5. Overlay the same points on the generalized-eigenfunction modulus.
%
%  REQUIRED PRECOMPUTATION AND LOCATION
%      Run ../duffing_regime_I.m, then duffing_upo_search.m.
%
%  RIGGEDDMD IMPLEMENTATION
%  The riggedDMD code used by the companion analysis is taken from:
%      https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
%  OUTPUTS
%      duffing_phase_1.pdf
%      duffing_mod_1.pdf
%      duffing_phase_2.pdf
%      duffing_mod_2.pdf
%      duffing_phase_3.pdf
%      duffing_mod_3.pdf
%
%  MAIN USER PARAMETERS
%      plotPeriods     minimal map periods to display
%      markerSize      size of the overlaid periodic points
%      saveFigures     save publication-quality PDF files
%
%  NOTES
%  The UPO search and plotting stages are intentionally separate: shooting
%  is comparatively expensive, whereas this script can quickly regenerate
%  figures after cosmetic changes.  The region labels are defined by the
%  phase convention saved by duffing_regime_I.m.  Legacy period-two records
%  containing adaptive ODE output are reduced to their two stroboscopic
%  endpoints before plotting.  This file is intended to reside in
%  duffing_utils/ beneath the repository root.
%
%  AUTHOR
%  Jason J. Bramburger
% -------------------------------------------------------------------------

clear; close all; clc;

%% User parameters and paths
plotPeriods = 1:3;
markerSize = 100;
saveFigures = true;

scriptDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(scriptDirectory);
outputDirectory = repositoryDirectory;
regimeFile = resolve_path(repositoryDirectory, ...
    'duffing_section_I_cache.mat');
upoFile = resolve_first_existing({ ...
    fullfile(repositoryDirectory,'duffing_upo_data.mat'), ...
    fullfile(repositoryDirectory,'duffing_results','duffing_upo_data.mat'), ...
    fullfile(scriptDirectory,'duffing_upo_data.mat')});

regimeData = load(regimeFile,'sectionData','region','logModulus');
upoData = load(upoFile,'orbitPeriods','orbitStrobePoints');

requiredRegimeFields = {'sectionData','region','logModulus'};
for fieldIndex = 1:numel(requiredRegimeFields)
    assert(isfield(regimeData,requiredRegimeFields{fieldIndex}), ...
        '%s does not contain %s.',regimeFile,requiredRegimeFields{fieldIndex})
end
assert(isfield(upoData,'orbitPeriods') && isfield(upoData,'orbitStrobePoints'), ...
    '%s does not contain the required orbit data.',upoFile)

sectionData = regimeData.sectionData;
region = regimeData.region(:);
logModulus = regimeData.logModulus(:);
orbitPeriods = upoData.orbitPeriods(:);
orbitStrobePoints = upoData.orbitStrobePoints;

assert(size(sectionData,1)==numel(region) && numel(region)==numel(logModulus), ...
    'The saved section, phase partition, and modulus have inconsistent lengths.')
assert(iscell(orbitStrobePoints) && ...
    numel(orbitStrobePoints)==numel(orbitPeriods), ...
    'The saved orbit periods and stroboscopic-point records are inconsistent.')

%% Publication style
regionColours = [27 158 119; 217 95 2; 117 112 179]/255;
orbitColour = [0 0 0];
pointSize = 10;

for mapPeriod = plotPeriods
    orbitIndices = find(orbitPeriods==mapPeriod);
    if isempty(orbitIndices)
        warning('No orbit of minimal period %d was found; skipping it.',mapPeriod)
        continue
    end
    [periodPoints,repairedRecordCount] = collect_period_points( ...
        orbitStrobePoints,orbitPeriods,orbitIndices);
    if repairedRecordCount > 0
        fprintf(['Reduced %d legacy period-two record(s) to their two ' ...
            'stroboscopic endpoints.\n'],repairedRecordCount)
    end
    fprintf('Minimal period %d: plotting %d section points from %d orbit(s).\n', ...
        mapPeriod,size(periodPoints,1),numel(orbitIndices))

    %% Phase partition with period-mapPeriod orbit intersections
    phaseFigure = figure('Name',sprintf('Duffing phase with period %d UPOs',mapPeriod), ...
        'Color','w','Position',[100 100 850 650]);
    phaseAxes = axes(phaseFigure);
    hold(phaseAxes,'on')
    phaseHandles = gobjects(3,1);
    for regionIndex = 1:3
        inRegion = region==regionIndex;
        phaseHandles(regionIndex) = scatter(phaseAxes,sectionData(inRegion,1), ...
            sectionData(inRegion,2),pointSize,regionColours(regionIndex,:), ...
            'filled','MarkerFaceAlpha',0.65,'MarkerEdgeAlpha',0.65);
    end
    orbitHandle = scatter(phaseAxes,periodPoints(:,1),periodPoints(:,2),markerSize, ...
        orbitColour,'filled','MarkerEdgeColor','w','LineWidth',1.2);
    xlabel(phaseAxes,'$x$','Interpreter','latex','FontSize',24)
    ylabel(phaseAxes,'$\dot{x}$','Interpreter','latex','FontSize',24)
    title(phaseAxes,sprintf('Minimal period $%d$',mapPeriod), ...
        'Interpreter','latex','FontSize',20)
    legend(phaseAxes,[phaseHandles; orbitHandle], ...
        {'$R_1$','$R_2$','$R_3$',sprintf('period-$%d$ UPO',mapPeriod)}, ...
        'Interpreter','latex','Location','best','FontSize',14)
    set(phaseAxes,'FontSize',16,'Layer','top')
    box(phaseAxes,'on')
    axis(phaseAxes,'tight')
    grid on

    %% Log-modulus with the same orbit intersections
    modulusFigure = figure('Name',sprintf('Duffing modulus with period %d UPOs',mapPeriod), ...
        'Color','w','Position',[100 100 850 650]);
    modulusAxes = axes(modulusFigure);
    scatter(modulusAxes,sectionData(:,1),sectionData(:,2),pointSize,logModulus, ...
        'filled','MarkerFaceAlpha',0.72,'MarkerEdgeAlpha',0.72)
    hold(modulusAxes,'on')
    scatter(modulusAxes,periodPoints(:,1),periodPoints(:,2),markerSize,orbitColour, ...
        'filled','MarkerEdgeColor','w','LineWidth',1.2)
    xlabel(modulusAxes,'$x$','Interpreter','latex','FontSize',24)
    ylabel(modulusAxes,'$\dot{x}$','Interpreter','latex','FontSize',24)
    title(modulusAxes,sprintf('Minimal period $%d$',mapPeriod), ...
        'Interpreter','latex','FontSize',20)
    colourBar = colorbar(modulusAxes);
    colourBar.Label.String = '$\log|\phi_{\varepsilon,\theta}|$';
    colourBar.Label.Interpreter = 'latex';
    colourBar.Label.FontSize = 20;
    colormap(modulusAxes,parula)
    set(modulusAxes,'FontSize',16,'Layer','top')
    box(modulusAxes,'on')
    axis(modulusAxes,'tight')
    grid on

    if saveFigures
        save_figure_pdf(phaseFigure,fullfile(outputDirectory, ...
            sprintf('duffing_phase_%d.pdf',mapPeriod)))
        save_figure_pdf(modulusFigure,fullfile(outputDirectory, ...
            sprintf('duffing_mod_%d.pdf',mapPeriod)))
    end
end

fprintf('Generated UPO overlays for minimal periods: %s\n', ...
    strjoin(string(intersect(plotPeriods,unique(orbitPeriods).')) ,', '));

%% Local functions
function [periodPoints,repairedRecordCount] = collect_period_points( ...
        orbitStrobePoints,orbitPeriods,orbitIndices)
%COLLECT_PERIOD_POINTS Validate and combine the requested section points.

pointCells = cell(numel(orbitIndices),1);
repairedRecordCount = 0;

for localIndex = 1:numel(orbitIndices)
    orbitIndex = orbitIndices(localIndex);
    expectedPointCount = orbitPeriods(orbitIndex);
    points = orbitStrobePoints{orbitIndex};

    assert(isnumeric(points) && ismatrix(points) && size(points,2)==2 && ...
        all(isfinite(points(:))), ...
        'Orbit %d does not contain a finite N-by-2 point array.',orbitIndex)

    if size(points,1)==expectedPointCount
        pointCells{localIndex} = points;
    elseif expectedPointCount==2 && size(points,1)>2
        % Older searches passed [0,T] to ode45.  For a two-entry time span,
        % ode45 returns its adaptive internal steps; only the endpoints are
        % successive stroboscopic intersections.
        pointCells{localIndex} = points([1 end],:);
        repairedRecordCount = repairedRecordCount+1;
    else
        error(['Orbit %d is tagged as minimal period %d but contains %d ' ...
            'stored points. Regenerate duffing_upo_data.mat with the ' ...
            'current duffing_upo_search.m.'], ...
            orbitIndex,expectedPointCount,size(points,1))
    end
end

periodPoints = vertcat(pointCells{:});
end

function save_figure_pdf(figHandle,fileName)
%SAVE_FIGURE_PDF Save a figure to a tightly sized vector PDF.

[fileDirectory,~,extension] = fileparts(fileName);
if isempty(extension)
    fileName = [fileName '.pdf'];
elseif ~strcmpi(extension,'.pdf')
    error('Figure output must have a .pdf extension.')
end

assert(isempty(fileDirectory) || isfolder(fileDirectory), ...
    'The figure output directory does not exist: %s',fileDirectory)

set(figHandle,'PaperUnits','centimeters')
set(figHandle,'Units','centimeters')
position = get(figHandle,'Position');
set(figHandle,'PaperSize',[position(3) position(4)])
set(figHandle,'PaperPositionMode','manual')
set(figHandle,'PaperPosition',[0 0 position(3) position(4)])
print(figHandle,fileName,'-dpdf','-painters')
end

function path = resolve_path(scriptDirectory,fileName)
candidates = {fullfile(scriptDirectory,fileName),fileName};
for candidateIndex = 1:numel(candidates)
    if isfile(candidates{candidateIndex})
        path = candidates{candidateIndex};
        return
    end
end
error('Could not find %s.',fileName);
end

function path = resolve_first_existing(candidates)
for candidateIndex = 1:numel(candidates)
    if isfile(candidates{candidateIndex})
        path = candidates{candidateIndex};
        return
    end
end
error('Could not find duffing_upo_data.mat.');
end
