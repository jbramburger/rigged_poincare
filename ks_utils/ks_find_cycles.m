% ============================================================
% DETECT CANDIDATE PERIODIC POINTS IN THE KURAMOTO–SIVASHINSKY
% POINCARE SECTION
% ============================================================
%
% DESCRIPTION:
% This script searches for candidate period-p points of the
% second-return map on a selected branch of the Kuramoto–Sivashinsky
% Poincaré section. It does so by computing the p-step return error
%
%       ||Ybranch(n+p,:) - Ybranch(n,:)||
%
% for each point along the branch and identifying those with the
% smallest return errors as candidate periodic points.
%
% Since several nearby time indices may correspond to the same
% geometric periodic orbit, the script then clusters the strongest
% candidates using a geometric tolerance and returns one
% representative point from each cluster.
%
% REQUIRED INPUTS (must exist in workspace):
%   Ybranch          : N x d array of Poincaré section points
%
%   NOTE:
%   For the Kuramoto–Sivashinsky application, Ybranch can be formed
%   from principal component projections (e.g., from
%   kuramoto_sivashinsky.m) as
%
%       Ybranch = [y1 y2 y3];
%
%   where y1, y2, and y3 are the leading principal component
%   coordinates of the Poincaré section data.
%
% USER-DEFINED PARAMETERS:
%   p                : target period for the second-return map
%   nCandidates      : number of lowest-return-error points retained
%   clusterTol       : geometric tolerance for clustering candidates
%
% OUTPUTS:
%   returnError      : p-step return error at each valid index
%   bestIndices      : strongest candidate indices before clustering
%   orbitReps        : representative indices after clustering
%
% VISUALIZATION:
%   - Scatter plot of the branch in principal component coordinates
%   - Candidate period-p points highlighted in red
%
% USAGE:
%   1. Ensure Ybranch is defined in the workspace.
%   2. For the KS case, construct Ybranch via
%
%          Ybranch = [y1 y2 y3];
%
%   3. Choose the target period p.
%   4. Adjust nCandidates and clusterTol if needed.
%   5. Run the script and inspect the highlighted candidates and
%      printed representative indices.
%
% REMARK:
% Small p-step return error indicates that a point is close to a
% period-p point of the second-return map. Because nearby samples
% along the same orbit may all have small return error, clustering
% is used to identify distinct geometric representatives.
% ============================================================

%% Find candidate period-p points in the KSE Poincare section

% --- User parameters ---
p = 5;                 % Target period for the second-return map
nCandidates = 50;      % Number of best return-error candidates to keep
clusterTol = 0.01;     % Geometric clustering tolerance

% --- Basic sizes ---
nPoints = size(Ybranch, 1);

% We can only evaluate the p-step return error up to index nPoints - p.
nValid = nPoints - p;

if nValid < 1
    error('Ybranch does not contain enough points for the requested period p = %d.', p);
end

% --- Step 1: compute p-step return errors ---
% returnError(n) = ||Ybranch(n+p,:) - Ybranch(n,:)||
returnError = NaN(nValid, 1);

for n = 1:nValid
    returnError(n) = norm(Ybranch(n + p, :) - Ybranch(n, :));
end

% --- Step 2: sort by smallest return error ---
[sortedErrors, sortedIndices] = sort(returnError, 'ascend');

% Keep only the strongest candidates
nKeep = min(nCandidates, numel(sortedIndices));
bestIndices = sortedIndices(1:nKeep);
bestErrors = sortedErrors(1:nKeep);

% --- Step 3: visualize the candidate points ---
figure(201)
clf
hold on
scatter(Ybranch(:,1), Ybranch(:,2), 8, [0.85 0.85 0.85], 'filled')
scatter(Ybranch(bestIndices,1), Ybranch(bestIndices,2), 50, 'r', 'filled')
xlabel('PC1')
ylabel('PC2')
title(sprintf('Candidate period-%d points', p))
axis equal
grid on
box on

% --- Step 4: cluster candidates geometrically ---
% Nearby candidate points may correspond to the same orbit. We keep
% one representative from each cluster.
used = false(nKeep, 1);
orbitReps = [];

for i = 1:nKeep
    if used(i)
        continue
    end

    % Start a new cluster with candidate i
    representativeIndex = bestIndices(i);
    representativePoint = Ybranch(representativeIndex, :);

    orbitReps(end+1, 1) = representativeIndex; %#ok<SAGROW>
    used(i) = true;

    % Mark all later candidates within clusterTol as belonging
    % to the same geometric cluster
    for j = i+1:nKeep
        candidatePoint = Ybranch(bestIndices(j), :);
        if norm(candidatePoint - representativePoint) < clusterTol
            used(j) = true;
        end
    end
end

% --- Step 5: report results ---
disp('Top candidate indices and p-step return errors:')
disp(table(bestIndices, bestErrors, ...
    'VariableNames', {'Index', 'ReturnError'}))

disp('Representative candidate indices after clustering:')
disp(orbitReps)

% Optional: also display the representative points themselves
disp('Representative candidate coordinates:')
disp(Ybranch(orbitReps, :))