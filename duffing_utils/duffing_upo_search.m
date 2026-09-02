% -------------------------------------------------------------------------
%  DUFFING_UPO_SEARCH  Search for unstable periodic orbits of the Regime I
%  Duffing stroboscopic map.
% -------------------------------------------------------------------------
%
%  DESCRIPTION
%  This script locates periodic points of the Poincare map using close
%  returns in a numerical trajectory as initial guesses, followed by a
%  shooting correction for the periodically forced Duffing equation.  It
%  searches all map periods from one through kMax, determines the minimal
%  period of each converged solution, removes duplicate copies (including
%  cyclic shifts along the same orbit), and saves the resulting orbit data.
%
%  WORKFLOW
%    1. Load stroboscopic data from duffing_psec_I.mat.
%    2. Identify the best k-step recurrences for k = 1,...,kMax.
%    3. Cluster nearby recurrences to obtain shooting initial guesses.
%    4. Solve F^k(z)-z = 0 with fsolve, or fminsearch when fsolve is absent.
%    5. Classify solutions by minimal period and deduplicate their cycles.
%    6. Save the orbit seeds and all section intersections for plotting.
%
%  MODEL
%      x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t),
%
%  with Regime I parameters
%      alpha = -1, beta = 1, delta = 0.3, omega = 1.2, gamma = 0.5.
%  The Poincare map samples the flow once per forcing period T=2*pi/omega.
%
%  RELATED RIGGEDDMD IMPLEMENTATION
%  The companion phase and modulus fields are computed using code from:
%      https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
%  INPUT
%      duffing_psec_I.mat containing x, an N-by-2 array whose columns are
%      position and velocity at successive forcing periods.
%
%  OUTPUT
%      ../duffing_results/duffing_upo_data.mat
%
%  MAIN USER PARAMETERS
%      kMax             largest map period to search
%      recurrenceFraction fraction of the closest k-step returns retained
%      clusterTolerance clustering radius in normalized section coordinates
%      shootingTolerance acceptance tolerance for F^k(z)-z
%      orbitDedupTolerance tolerance for identifying cyclic copies
%
%  NOTES
%  The Optimization Toolbox is optional.  When fsolve is unavailable, the
%  script minimizes the squared shooting residual with fminsearch.  Orbit
%  discovery can be computationally expensive; run duffing_upo_plot.m to
%  regenerate figures without repeating this search.  This file is intended
%  to reside in duffing_utils/ beneath the repository root.
%
%  AUTHOR
%  Jason J. Bramburger
% -------------------------------------------------------------------------

clear; close all; clc;

%% Paths and input data
scriptDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(scriptDirectory);

dataFile = resolve_path(scriptDirectory,'duffing_psec_I.mat');
data = load(dataFile,'x');
if ~isfield(data,'x') || size(data.x,2) ~= 2
    error('duffing_psec_I.mat must contain an N-by-2 variable named x.');
end
sectionData = data.x;

resultsDirectory = fullfile(repositoryDirectory,'duffing_results');
if ~exist(resultsDirectory,'dir')
    mkdir(resultsDirectory);
end

%% Duffing parameters and numerical tolerances
alpha = -1;
beta  = 1;
delta = 0.3;
omega = 1.2;
gamma = 0.5;
forcingPeriod = 2*pi/omega;

kMax                  = 3;
recurrenceFraction    = 0.005;
clusterTolerance      = 0.04;
shootingTolerance     = 1e-10;
minimalPeriodTolerance = 1e-7;
orbitDedupTolerance   = 1e-5;

odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-12);
haveFsolve = exist('fsolve','file') == 2;
if haveFsolve
    solveOptions = optimoptions('fsolve','Display','off', ...
        'FunctionTolerance',1e-14,'StepTolerance',1e-14, ...
        'OptimalityTolerance',1e-14,'MaxIterations',1000);
else
    solveOptions = optimset('Display','off','TolX',1e-13,'TolFun',1e-13, ...
        'MaxIter',2e4,'MaxFunEvals',2e4);
end

%% Normalize the section for recurrence detection
sectionMean = mean(sectionData,1);
sectionScale = std(sectionData,0,1);
sectionScale(sectionScale == 0) = 1;
normalizedSection = (sectionData-sectionMean)./sectionScale;
numberOfPoints = size(sectionData,1);

candidateSeeds = zeros(0,2);
candidatePeriods = zeros(0,1);
candidateResiduals = zeros(0,1);

%% Search periods 1,...,kMax
for mapPeriod = 1:kMax
    fprintf('\nSearching for period-%d points\n',mapPeriod);
    if numberOfPoints <= mapPeriod
        continue
    end

    recurrenceDistance = vecnorm( ...
        normalizedSection(1:end-mapPeriod,:)- ...
        normalizedSection(1+mapPeriod:end,:),2,2);
    cutoff = quantile(recurrenceDistance,recurrenceFraction);
    recurrenceIndices = find(recurrenceDistance <= cutoff);
    clusterCenters = greedy_cluster( ...
        normalizedSection(recurrenceIndices,:),clusterTolerance);
    physicalCenters = clusterCenters.*sectionScale+sectionMean;
    fprintf('  recurrence clusters: %d\n',size(physicalCenters,1));

    for centerIndex = 1:size(physicalCenters,1)
        initialGuess = physicalCenters(centerIndex,:).';
        shootingTime = mapPeriod*forcingPeriod;
        residualFunction = @(z) shooting_residual(z,shootingTime, ...
            alpha,beta,delta,omega,gamma,odeOptions);

        if haveFsolve
            correctedSeed = fsolve(residualFunction,initialGuess,solveOptions);
        else
            objective = @(z) norm(residualFunction(z(:)))^2;
            correctedSeed = fminsearch(objective,initialGuess,solveOptions);
        end
        correctedSeed = correctedSeed(:);
        residualNorm = norm(residualFunction(correctedSeed));
        if ~all(isfinite(correctedSeed)) || residualNorm >= shootingTolerance
            continue
        end

        minimalPeriod = mapPeriod;
        for testPeriod = 1:mapPeriod-1
            testResidual = shooting_residual(correctedSeed, ...
                testPeriod*forcingPeriod,alpha,beta,delta,omega,gamma,odeOptions);
            if norm(testResidual) < minimalPeriodTolerance
                minimalPeriod = testPeriod;
                residualNorm = norm(testResidual);
                break
            end
        end

        candidateSeeds(end+1,:) = correctedSeed.'; %#ok<SAGROW>
        candidatePeriods(end+1,1) = minimalPeriod; %#ok<SAGROW>
        candidateResiduals(end+1,1) = residualNorm; %#ok<SAGROW>
    end
end

if isempty(candidateSeeds)
    error(['No periodic orbits converged. Increase recurrenceFraction or ', ...
        'clusterTolerance, or relax shootingTolerance.']);
end

%% Deduplicate complete cycles, including cyclic shifts
[candidateResiduals,order] = sort(candidateResiduals);
candidateSeeds = candidateSeeds(order,:);
candidatePeriods = candidatePeriods(order);

orbitSeeds = zeros(0,2);
orbitPeriods = zeros(0,1);
orbitResiduals = zeros(0,1);
orbitStrobePoints = cell(0,1);

for candidateIndex = 1:size(candidateSeeds,1)
    seed = candidateSeeds(candidateIndex,:).';
    period = candidatePeriods(candidateIndex);
    strobePoints = generate_strobe_points(seed,period,forcingPeriod, ...
        alpha,beta,delta,omega,gamma,odeOptions);
    assert(size(strobePoints,1)==period && size(strobePoints,2)==2, ...
        'A period-%d orbit did not produce exactly %d section points.', ...
        period,period)

    duplicate = false;
    for orbitIndex = 1:numel(orbitPeriods)
        if orbitPeriods(orbitIndex) ~= period
            continue
        end
        pairwiseDistance = vecnorm(orbitStrobePoints{orbitIndex}-seed.',2,2);
        if min(pairwiseDistance) < orbitDedupTolerance
            duplicate = true;
            break
        end
    end
    if duplicate
        continue
    end

    orbitSeeds(end+1,:) = seed.'; %#ok<SAGROW>
    orbitPeriods(end+1,1) = period; %#ok<SAGROW>
    orbitResiduals(end+1,1) = candidateResiduals(candidateIndex); %#ok<SAGROW>
    orbitStrobePoints{end+1,1} = strobePoints; %#ok<SAGROW>
end

allStrobePoints = vertcat(orbitStrobePoints{:});
strobePeriodCells = arrayfun(@(j) ...
    repmat(orbitPeriods(j),orbitPeriods(j),1), ...
    (1:numel(orbitPeriods)).','UniformOutput',false);
allStrobePeriods = vertcat(strobePeriodCells{:});

orbitNumber = (1:numel(orbitPeriods)).';
summaryTable = table(orbitNumber,orbitPeriods,orbitResiduals, ...
    'VariableNames',{'Orbit','MinimalPeriod','ShootingResidual'});
disp(summaryTable);

outputFile = fullfile(resultsDirectory,'duffing_upo_data.mat');
save(outputFile,'orbitSeeds','orbitPeriods','orbitResiduals', ...
    'orbitStrobePoints','allStrobePoints','allStrobePeriods', ...
    'summaryTable','alpha','beta','delta','omega','gamma', ...
    'forcingPeriod','kMax','recurrenceFraction','clusterTolerance', ...
    'shootingTolerance','minimalPeriodTolerance','orbitDedupTolerance');
fprintf('\nSaved %d distinct orbit(s) to %s\n',numel(orbitPeriods),outputFile);

%% Local functions
function centers = greedy_cluster(points,tolerance)
centers = zeros(0,size(points,2));
members = cell(0,1);
for pointIndex = 1:size(points,1)
    point = points(pointIndex,:);
    if isempty(centers)
        centers = point;
        members{1,1} = point;
        continue
    end
    [minimumDistance,clusterIndex] = min(vecnorm(centers-point,2,2));
    if minimumDistance < tolerance
        members{clusterIndex}(end+1,:) = point;
        centers(clusterIndex,:) = mean(members{clusterIndex},1);
    else
        centers(end+1,:) = point; %#ok<AGROW>
        members{end+1,1} = point; %#ok<AGROW>
    end
end
end

function residual = shooting_residual(seed,integrationTime,alpha,beta,delta,omega,gamma,options)
[~,trajectory] = ode45(@(t,z) duffing_rhs(t,z,alpha,beta,delta,omega,gamma), ...
    [0 integrationTime],seed(:),options);
residual = trajectory(end,:).'-seed(:);
end

function points = generate_strobe_points(seed,period,forcingPeriod,alpha,beta,delta,omega,gamma,options)
% Return exactly one point at each forcing period.  In particular, avoid
% passing [0,T] as the ode45 output vector for period two: MATLAB interprets
% a two-entry vector as an integration interval and returns adaptive internal
% steps rather than only the requested stroboscopic endpoints.
points = zeros(period,numel(seed));
points(1,:) = seed(:).';
state = seed(:);
currentTime = 0;

for pointIndex = 2:period
    nextTime = currentTime+forcingPeriod;
    [~,trajectory] = ode45( ...
        @(t,z) duffing_rhs(t,z,alpha,beta,delta,omega,gamma), ...
        [currentTime nextTime],state,options);
    state = trajectory(end,:).';
    points(pointIndex,:) = state.';
    currentTime = nextTime;
end
end

function dz = duffing_rhs(t,z,alpha,beta,delta,omega,gamma)
dz = [z(2); gamma*cos(omega*t)-delta*z(2)-alpha*z(1)-beta*z(1)^3];
end

function path = resolve_path(scriptDirectory,fileName)
candidates = {fullfile(scriptDirectory,fileName), ...
    fullfile(scriptDirectory,'data',fileName),fileName};
for candidateIndex = 1:numel(candidates)
    if isfile(candidates{candidateIndex})
        path = candidates{candidateIndex};
        return
    end
end
error('Could not find %s.',fileName);
end