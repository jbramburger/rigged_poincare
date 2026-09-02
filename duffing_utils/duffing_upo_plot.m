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
%  phase convention saved by duffing_regime_I.m.  This file is intended to
%  reside in duffing_utils/ beneath the repository root.
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

%% Publication style
regionColours = [36 122 254; 230 159 0; 0 120 0]/255;
orbitColour = [0 0 0];
pointSize = 10;

for mapPeriod = plotPeriods
    orbitIndices = find(orbitPeriods==mapPeriod);
    if isempty(orbitIndices)
        warning('No orbit of minimal period %d was found; skipping it.',mapPeriod)
        continue
    end
    periodPoints = vertcat(orbitStrobePoints{orbitIndices});

    %% Phase partition with period-mapPeriod orbit intersections
    phaseFigure = figure('Name',sprintf('Duffing phase with period %d UPOs',mapPeriod), ...
        'Color','w','Position',[100 100 850 650]);
    hold on
    phaseHandles = gobjects(3,1);
    for regionIndex = 1:3
        inRegion = region==regionIndex;
        phaseHandles(regionIndex) = scatter(sectionData(inRegion,1), ...
            sectionData(inRegion,2),pointSize,regionColours(regionIndex,:), ...
            'filled','MarkerFaceAlpha',0.65,'MarkerEdgeAlpha',0.65);
    end
    orbitHandle = scatter(periodPoints(:,1),periodPoints(:,2),markerSize, ...
        orbitColour,'filled','MarkerEdgeColor','w','LineWidth',1.2);
    xlabel('$x$','Interpreter','latex','FontSize',24)
    ylabel('$\dot{x}$','Interpreter','latex','FontSize',24)
    title(sprintf('Minimal period $%d$',mapPeriod), ...
        'Interpreter','latex','FontSize',20)
    legend([phaseHandles; orbitHandle], ...
        {'$R_1$','$R_2$','$R_3$',sprintf('period-$%d$ UPO',mapPeriod)}, ...
        'Interpreter','latex','Location','best','FontSize',14)
    set(gca,'FontSize',16,'Layer','top')
    box on
    axis tight

    %% Log-modulus with the same orbit intersections
    modulusFigure = figure('Name',sprintf('Duffing modulus with period %d UPOs',mapPeriod), ...
        'Color','w','Position',[100 100 850 650]);
    scatter(sectionData(:,1),sectionData(:,2),pointSize,logModulus, ...
        'filled','MarkerFaceAlpha',0.72,'MarkerEdgeAlpha',0.72)
    hold on
    scatter(periodPoints(:,1),periodPoints(:,2),markerSize,orbitColour, ...
        'filled','MarkerEdgeColor','w','LineWidth',1.2)
    xlabel('$x$','Interpreter','latex','FontSize',24)
    ylabel('$\dot{x}$','Interpreter','latex','FontSize',24)
    title(sprintf('Minimal period $%d$',mapPeriod), ...
        'Interpreter','latex','FontSize',20)
    colourBar = colorbar;
    colourBar.Label.String = '$\log|\phi|$';
    colourBar.Label.Interpreter = 'latex';
    colourBar.Label.FontSize = 20;
    colormap(parula)
    set(gca,'FontSize',16,'Layer','top')
    box on
    axis tight

    if saveFigures
        exportgraphics(phaseFigure,fullfile(outputDirectory, ...
            sprintf('duffing_phase_%d.pdf',mapPeriod)),'ContentType','vector');
        exportgraphics(modulusFigure,fullfile(outputDirectory, ...
            sprintf('duffing_mod_%d.pdf',mapPeriod)),'ContentType','vector');
    end
end

fprintf('Generated UPO overlays for minimal periods: %s\n', ...
    strjoin(string(intersect(plotPeriods,unique(orbitPeriods).')) ,', '));

%% Local function
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
