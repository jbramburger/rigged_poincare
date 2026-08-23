% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE PERIODICALLY FORCED VAN DER POL POINCARE MAP
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This supplementary script analyzes the stroboscopic Poincare map of a
% periodically forced Van der Pol oscillator.  It applies rigged Dynamic
% Mode Decomposition (riggedDMD) to a delay-coordinate dictionary and uses
% the resulting generalized eigenfunction to examine phase transport and
% low-modulus defect regions on the section.
%
% The script computes:
%   1. The stroboscopic Poincare section
%   2. The riggedDMD spectral density
%   3. A smoothed generalized eigenfunction near theta = 0.81*pi
%   4. The phase return map and wrapped phase-increment distribution
%   5. Two empirical phase-transport branches
%   6. The relationship between branch membership and eigenfunction modulus
%   7. Return and transition statistics for a low-modulus defect set
%   8. Spectral robustness across x, xdot, and x+xdot observables
%
% MODEL
% The periodically forced Van der Pol equation is
%
%   x'' + mu*(x^2-1)*x' + x = forcingAmplitude*cos(forcingFrequency*t),
%
% with
%   mu = 3, forcingAmplitude = 5, forcingFrequency = 1.788.
%
% The Poincare map is obtained by sampling once per forcing period
% T = 2*pi/forcingFrequency.
%
% RIGGEDDMD IMPLEMENTATION
% The riggedDMD code used here is taken from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% INPUT
%   vdp_psec.mat              : optional cached Poincare data containing an
%                               N-by-2 array x = [position,velocity]
%   main_routines/            : riggedDMD and its dependencies
%
% OUTPUTS
%   vanderpol_results/vanderpol_results.mat
%   vanderpol_results/vanderpol_poincare_section.pdf
%   vanderpol_results/vanderpol_spectral_density.pdf
%   vanderpol_results/vanderpol_observable_robustness.pdf
%   vanderpol_results/vanderpol_eigenfunction.pdf
%   vanderpol_results/vanderpol_phase_transport.pdf
%   vanderpol_results/vanderpol_transport_branches.pdf
%   vanderpol_results/vanderpol_defect_set.pdf
%
% MAIN USER PARAMETERS
%   regenerateData            : regenerate and overwrite vdp_psec.mat
%   nTransient                : discarded stroboscopic samples
%   nData                     : retained stroboscopic samples
%   dim                       : number of delays
%   epsilon                   : riggedDMD smoothing parameter
%   thetaTarget               : generalized-eigenfunction frequency
%   branchThreshold           : phase-increment branch threshold in radians
%   defectThreshold           : threshold applied to log|phi|
%   saveFigures               : save PDF figures
%
% NOTES
%   - This example records supplementary numerical tests and is not part of
%     the accompanying paper.
%   - If vdp_psec.mat is absent, the script generates and caches the section
%     data automatically.
%   - The branch and defect thresholds reproduce the exploratory tests in
%     the original script and are exposed as user parameters.
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

dataFile = 'vdp_psec.mat';
routinesDirectory = 'main_routines';
outputDirectory = 'vanderpol_results';

regenerateData = false;
nTransient = 200;
nData = 20000;
initialCondition = [1;0];

mu = 3;
forcingAmplitude = 5;
forcingFrequency = 1.788;

dim = 20;
epsilon = 0.25;
order = 2;
thetaTarget = 0.81*pi;
thetaGrid = linspace(-pi,pi,1601);

branchThreshold = 1.5;
defectThreshold = -1;
nestedDefectThresholds = [-1 -2 -3];
runObservableRobustness = true;
saveFigures = true;

odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-12);

%% Locate dependencies and load or generate the Poincare data

scriptDirectory = fileparts(mfilename('fullpath'));
if isempty(scriptDirectory)
    scriptDirectory = pwd;
end

routinesPath = resolve_path(scriptDirectory,routinesDirectory,'directory');
addpath(routinesPath)
assert(exist('riggedDMD','file') == 2, ...
    'riggedDMD.m was not found in %s.',routinesPath)

dataPath = fullfile(scriptDirectory,dataFile);
if regenerateData || ~isfile(dataPath)
    fprintf('Generating %d forced Van der Pol Poincare points...\n',nData)
    x = generate_poincare_data(mu,forcingAmplitude,forcingFrequency, ...
        initialCondition,nTransient,nData,odeOptions);
    save(dataPath,'x','mu','forcingAmplitude','forcingFrequency', ...
        'nTransient','nData','initialCondition')
else
    S = load(dataPath);
    assert(isfield(S,'x'),'vdp_psec.mat must contain the variable x.')
    x = S.x;
end

assert(isnumeric(x) && size(x,2) == 2, ...
    'x must be an N-by-2 array [position, velocity].')
assert(size(x,1) > dim+1,'The Poincare dataset is too short for dim delays.')
assert(all(isfinite(x(:))),'The Poincare data contain NaN or Inf values.')

%% Construct the delay-coordinate dictionary

observable = normalize_observable(x(:,1));
[PsiX,PsiY] = delay_dictionary(observable,dim);

nSnapshots = size(PsiX,1);
sectionData = x(1:nSnapshots,:);
quadratureWeight = 1/size(x,1);

gCoeffs = zeros(dim,1);
gCoeffs(1) = 1;

%% Apply riggedDMD

[gMode,spectralDensity] = run_rigged_dmd(PsiX,PsiY,quadratureWeight, ...
    epsilon,thetaTarget,order,gCoeffs,thetaGrid);

phi = PsiX*gMode;
phase = angle(phi);
modulus = abs(phi);
logModulus = log(max(modulus,realmin));

[peakLocation,peakHeight] = local_maximum( ...
    thetaGrid,real(spectralDensity),thetaTarget,0.10*pi);

%% Phase transport and low-modulus defect diagnostics

phaseCurrent = phase(1:end-1);
phaseNext = phase(2:end);
phaseIncrement = wrap_to_pi(phaseNext-phaseCurrent);

transportBranch = phaseIncrement > branchThreshold;
branchFraction = [mean(~transportBranch) mean(transportBranch)];
stepLogModulus = logModulus(1:end-1);
branchMeanLogModulus = [mean(stepLogModulus(~transportBranch)) ...
    mean(stepLogModulus(transportBranch))];

defect = logModulus < defectThreshold;
defectCurrent = defect(1:end-1);
defectNext = defect(2:end);
defectFraction = mean(defect);

defectTransitionCounts = zeros(2,2);
defectTransitionCounts(1,1) = sum(~defectCurrent & ~defectNext);
defectTransitionCounts(1,2) = sum(~defectCurrent &  defectNext);
defectTransitionCounts(2,1) = sum( defectCurrent & ~defectNext);
defectTransitionCounts(2,2) = sum( defectCurrent &  defectNext);
defectTransitionMatrix = defectTransitionCounts ./ ...
    max(sum(defectTransitionCounts,2),1);

phaseStdOnDefect = std(phaseIncrement(defectCurrent));
phaseStdOffDefect = std(phaseIncrement(~defectCurrent));
defectIndices = find(defect);
defectReturnGaps = diff(defectIndices);

fprintf('\nForced Van der Pol riggedDMD diagnostics\n')
fprintf('-----------------------------------------\n')
fprintf('Poincare points:                    %d\n',size(x,1))
fprintf('Delay dimension:                    %d\n',dim)
fprintf('Smoothing parameter:                %.4f\n',epsilon)
fprintf('Target theta/pi:                    %.5f\n',thetaTarget/pi)
fprintf('Local spectral maximum theta/pi:    %.5f\n',peakLocation/pi)
fprintf('Local spectral maximum height:      %.4e\n',peakHeight)
fprintf('Branch fractions below/above cut:   %.4f / %.4f\n',branchFraction)
fprintf('Mean log modulus by branch:         %.4f / %.4f\n',branchMeanLogModulus)
fprintf('Defect threshold:                   %.4f\n',defectThreshold)
fprintf('Defect fraction:                    %.4f\n',defectFraction)
fprintf('std(phase increment) on defect:     %.4f\n',phaseStdOnDefect)
fprintf('std(phase increment) off defect:    %.4f\n',phaseStdOffDefect)
disp('Row-normalized defect transition matrix:')
disp(defectTransitionMatrix)

%% Spectral robustness across observables

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
            local_maximum(thetaGrid,observableSpectra(:,j),thetaTarget,0.10*pi);
    end

    fprintf('\nObservable robustness near theta/pi = %.3f\n',thetaTarget/pi)
    for j = 1:numel(observableNames)
        fprintf('  %-8s: peak theta/pi = %.5f, height = %.4e\n', ...
            observableNames{j},observablePeakLocation(j)/pi,observablePeakHeight(j))
    end
end

%% Plot the Poincare section

figure('Name','Van der Pol Poincare section','Color','w')
scatter(x(:,1),x(:,2),10,'k','filled','MarkerFaceAlpha',0.65)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\dot{x}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
axis tight
grid on
box on

%% Plot the spectral density

figure('Name','Van der Pol spectral density','Color','w')
plot(thetaGrid/pi,real(spectralDensity),'Color',[1 69/255 79/255], ...
    'LineWidth',3)
xline(thetaTarget/pi,'k--','target','LineWidth',1.5)
xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('Spectral density','Interpreter','latex','FontSize',20)
set(gca,'FontSize',16)
xlim([0 1])
grid on
box on

%% Plot observable robustness

if runObservableRobustness
    figure('Name','Van der Pol observable robustness','Color','w')
    hold on
    colours = lines(numel(observableNames));
    for j = 1:numel(observableNames)
        plot(thetaGrid/pi,observableSpectra(:,j),'LineWidth',2, ...
            'Color',colours(j,:),'DisplayName',observableNames{j})
    end
    xline(thetaTarget/pi,'k--','target','LineWidth',1.5, ...
        'HandleVisibility','off')
    xlabel('$\theta/\pi$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    ylabel('Spectral density','Interpreter','latex','FontSize',20)
    set(gca,'FontSize',16)
    xlim([0 1])
    legend('Location','best')
    grid on
    box on
end

%% Plot the generalized eigenfunction

figure('Name','Van der Pol generalized eigenfunction','Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
scatter(sectionData(:,1),sectionData(:,2),10,phase,'filled')
xlabel('$x$','Interpreter','latex','FontSize',20)
ylabel('$\dot{x}$','Interpreter','latex','FontSize',20)
title('Phase','FontSize',18)
set(gca,'FontSize',14)
axis tight
box on
colorbar

nexttile
scatter(sectionData(:,1),sectionData(:,2),10,logModulus,'filled')
xlabel('$x$','Interpreter','latex','FontSize',20)
ylabel('$\dot{x}$','Interpreter','latex','FontSize',20)
title('$\log|\phi|$','Interpreter','latex','FontSize',18)
set(gca,'FontSize',14)
axis tight
box on
colorbar
colormap(parula)

%% Plot phase transport

figure('Name','Van der Pol phase transport','Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
scatter(phaseCurrent,phaseNext,8,phaseIncrement,'filled')
xlabel('$\psi_n$','Interpreter','latex','FontSize',20)
ylabel('$\psi_{n+1}$','Interpreter','latex','FontSize',20)
title('Phase return map','FontSize',18)
set(gca,'FontSize',14)
grid on
box on
colorbar

nexttile
histogram(phaseIncrement,100,'FaceColor',[36 122 254]/255, ...
    'EdgeColor','none')
xline(branchThreshold,'k--','branch threshold','LineWidth',1.5)
xlabel('$\Delta\psi_n$','Interpreter','latex','FontSize',20)
ylabel('Count','FontSize',18)
title('Wrapped phase increment','FontSize',18)
set(gca,'FontSize',14)
grid on
box on

%% Plot transport branches on current and next section points

figure('Name','Van der Pol transport branches','Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
scatter(sectionData(1:end-1,1),sectionData(1:end-1,2),10, ...
    double(transportBranch),'filled')
xlabel('$x_n$','Interpreter','latex','FontSize',20)
ylabel('$\dot{x}_n$','Interpreter','latex','FontSize',20)
title('Current point','FontSize',18)
set(gca,'FontSize',14)
axis tight
box on

nexttile
scatter(sectionData(2:end,1),sectionData(2:end,2),10, ...
    double(transportBranch),'filled')
xlabel('$x_{n+1}$','Interpreter','latex','FontSize',20)
ylabel('$\dot{x}_{n+1}$','Interpreter','latex','FontSize',20)
title('Next point','FontSize',18)
set(gca,'FontSize',14)
axis tight
box on
colormap(lines(2))

%% Plot the low-modulus defect set and its return gaps

figure('Name','Van der Pol defect set','Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
scatter(sectionData(:,1),sectionData(:,2),7,[0.82 0.82 0.82],'filled')
hold on
defectColours = [36 122 254; 230 159 0; 0 0 0]/255;
defectSizes = [12 18 24];
for j = 1:numel(nestedDefectThresholds)
    inNestedDefect = logModulus < nestedDefectThresholds(j);
    scatter(sectionData(inNestedDefect,1),sectionData(inNestedDefect,2), ...
        defectSizes(j),defectColours(j,:),'filled')
end
xlabel('$x$','Interpreter','latex','FontSize',20)
ylabel('$\dot{x}$','Interpreter','latex','FontSize',20)
title('Nested low-modulus regions','FontSize',18)
legend({'all points','$\log|\phi|<-1$','$\log|\phi|<-2$', ...
    '$\log|\phi|<-3$'},'Interpreter','latex','Location','best')
set(gca,'FontSize',14)
axis tight
box on

nexttile
if isempty(defectReturnGaps)
    text(0.5,0.5,'Fewer than two defect visits','HorizontalAlignment','center')
    axis off
else
    histogram(defectReturnGaps,50,'FaceColor',[230 159 0]/255, ...
        'EdgeColor','none')
    xlabel('Return gap','FontSize',18)
    ylabel('Count','FontSize',18)
    title('Returns to the defect set','FontSize',18)
    set(gca,'FontSize',14)
    grid on
    box on
end

%% Save results and figures

outputPath = fullfile(scriptDirectory,outputDirectory);
if ~isfolder(outputPath)
    mkdir(outputPath)
end

resultsFile = fullfile(outputPath,'vanderpol_results.mat');
save(resultsFile,'x','sectionData','mu','forcingAmplitude', ...
    'forcingFrequency','dim','epsilon','order','thetaTarget','thetaGrid', ...
    'spectralDensity','peakLocation','peakHeight','phi','phase','modulus', ...
    'logModulus','phaseIncrement','branchThreshold','transportBranch', ...
    'branchFraction','branchMeanLogModulus','defectThreshold','defect', ...
    'defectFraction','defectTransitionCounts','defectTransitionMatrix', ...
    'phaseStdOnDefect','phaseStdOffDefect','defectReturnGaps', ...
    'observableNames','observablePeakLocation','observablePeakHeight')

if saveFigures
    save_named_figure('Van der Pol Poincare section',outputPath, ...
        'vanderpol_poincare_section')
    save_named_figure('Van der Pol spectral density',outputPath, ...
        'vanderpol_spectral_density')
    if runObservableRobustness
        save_named_figure('Van der Pol observable robustness',outputPath, ...
            'vanderpol_observable_robustness')
    end
    save_named_figure('Van der Pol generalized eigenfunction',outputPath, ...
        'vanderpol_eigenfunction')
    save_named_figure('Van der Pol phase transport',outputPath, ...
        'vanderpol_phase_transport')
    save_named_figure('Van der Pol transport branches',outputPath, ...
        'vanderpol_transport_branches')
    save_named_figure('Van der Pol defect set',outputPath, ...
        'vanderpol_defect_set')
end

fprintf('\nSaved forced Van der Pol results to %s\n',resultsFile)

%% Local functions

function x = generate_poincare_data(mu,forcingAmplitude,forcingFrequency, ...
        initialCondition,nTransient,nData,odeOptions)
    forcingPeriod = 2*pi/forcingFrequency;
    sampleTimes = (0:nTransient+nData)*forcingPeriod;
    [~,trajectory] = ode45(@(t,z) forced_vdp_rhs(t,z,mu, ...
        forcingAmplitude,forcingFrequency),sampleTimes,initialCondition,odeOptions);
    x = trajectory(nTransient+2:end,:);
end

function dz = forced_vdp_rhs(t,z,mu,forcingAmplitude,forcingFrequency)
    dz = [z(2); mu*(1-z(1)^2)*z(2)-z(1)+ ...
        forcingAmplitude*cos(forcingFrequency*t)];
end

function g = normalize_observable(g)
    g = g(:);
    scale = std(g);
    assert(scale > 0,'The selected observable has zero variance.')
    g = (g-mean(g))/scale;
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

function [location,height] = local_maximum(thetaGrid,density,target,halfWidth)
    window = abs(wrap_to_pi(thetaGrid-target)) <= halfWidth;
    indices = find(window);
    [height,localIndex] = max(real(density(window)));
    location = thetaGrid(indices(localIndex));
end

function pathOut = resolve_path(scriptDirectory,name,pathType)
    pathOut = fullfile(scriptDirectory,name);
    if strcmp(pathType,'file') && isfile(pathOut)
        return
    elseif strcmp(pathType,'directory') && isfolder(pathOut)
        return
    end
    if strcmp(pathType,'file') && isfile(name)
        pathOut = name;
        return
    elseif strcmp(pathType,'directory') && isfolder(name)
        pathOut = name;
        return
    end
    error('Could not locate %s.',name)
end

function save_named_figure(figureName,outputPath,fileName)
    figHandle = findobj(groot,'Type','figure','Name',figureName);
    assert(~isempty(figHandle),'Could not find figure named %s.',figureName)
    save_figure_pdf(figHandle(1),fullfile(outputPath,[fileName '.pdf']))
end

function save_figure_pdf(figHandle,fileName)
    try
        exportgraphics(figHandle,fileName,'ContentType','vector')
    catch
        set(figHandle,'PaperPositionMode','auto')
        print(figHandle,fileName,'-dpdf','-painters')
    end
end
