% ============================================================
% MODULUS-BASED DETECTION OF PERIOD-3 ORBITS AND PARTITIONING
% ============================================================
%
% DESCRIPTION:
% This script implements a robust, data-driven method for detecting
% a period-3 orbit and its preimages from a Poincaré section using
% Koopman eigenfunction data. The approach leverages local minima
% of the (smoothed) logarithmic modulus of an eigenfunction to
% identify dynamically significant points.
%
% The script proceeds in two main stages:
%
%   (1) Orbit Detection:
%       - Computes the smoothed log-modulus of a Koopman eigenfunction
%         along an ordered branch of the Poincaré section.
%       - Identifies well-separated local minima as candidate points.
%       - Filters candidates using a 3-step return test to isolate
%         a period-3 orbit.
%       - Detects preimages as points that map into the orbit in one
%         step but are not themselves periodic.
%
%   (2) Symbolic Partition and Transition Analysis:
%       - Uses the detected period-3 orbit to partition the ordered
%         branch into three regions.
%       - Assigns symbolic labels to each point in both ordered and
%         time coordinates.
%       - Constructs a 3×3 transition matrix approximating the
%         second-return map dynamics between regions.
%
% REQUIRED INPUTS (must exist in workspace):
%   C                : Koopman eigenfunction evaluated along trajectory
%   Ybranch          : N x d array of Poincaré section points
%   ord              : ordering of points along the branch
%
%   NOTE:
%   For the Kuramoto–Sivashinsky application, Ybranch can be formed
%   from principal component projections (e.g., from kuramoto_sivashinsky.m) as:
%
%       Ybranch = [y1 y2 y3];
%
%   where y1, y2, y3 are the leading principal component coordinates
%   of the Poincaré section data.
%
% KEY PARAMETERS (user adjustable):
%   w                : smoothing window for modulus
%   minSep           : minimum separation between detected minima
%   nKeep            : number of candidate minima retained
%   orbitTol         : tolerance for identifying period-3 points
%   preTol           : tolerance for identifying preimages
%
% OUTPUTS:
%   orbit_cycle      : indices of detected period-3 orbit
%   preimage_ids     : indices of candidate preimages
%   P3               : 3×3 transition matrix between regions
%   labels_ord_3     : region labels along ordered branch
%   labels_time_3    : region labels in time order
%
% VISUALIZATION:
%   - Smoothed modulus with detected minima (orbit in red, preimages in blue)
%   - Phase-space plot of detected orbit and preimages
%   - Partitioned branch colored by symbolic region
%   - Symbolic itinerary of the trajectory
%
% USAGE:
%   1. Ensure C, Ybranch, and ord are defined in the workspace.
%   2. (KS case) Construct Ybranch from PCA data via Ybranch = [y1 y2 y3].
%   3. Adjust parameters (w, minSep, orbitTol, etc.) if needed.
%   4. Run the script.
%   5. Inspect printed diagnostics, plots, and transition matrix P3.
%
% REMARK:
% This method is particularly effective for uncovering low-period
% structures embedded in high-dimensional Poincaré sections where
% the dynamics exhibit approximate low-dimensional organization.
% ============================================================
close all
clc

% PARAMETERS
w = 11;          % smoothing window for modulus
minSep = 20;     % minimum separation between minima along ordered branch
nKeep = 20;      % keep this many deepest minima
orbitTol = 0.1; % tolerance for 3-step return error (adjust)
preTol = 0.1;   % tolerance for 1-step image into orbit (adjust)

% Ordered smoothed modulus
modC = log(abs(C(:)));
modC_ord = modC(ord);
modC_smooth = movmean(modC_ord, w);

% Local minima detection
isMin = false(size(modC_smooth));
N = length(modC_smooth);
for i = 2:N-1
    if modC_smooth(i) <= modC_smooth(i-1) && modC_smooth(i) < modC_smooth(i+1)
        isMin(i) = true;
    elseif modC_smooth(i) < modC_smooth(i-1) && modC_smooth(i) <= modC_smooth(i+1)
        isMin(i) = true;
    end
end

min_idx_ord_all = find(isMin);

% Greedy separation, deepest first
[~, sortAll] = sort(modC_smooth(min_idx_ord_all), 'ascend');
selected = [];
for ii = 1:length(sortAll)
    idx0 = min_idx_ord_all(sortAll(ii));
    if isempty(selected) || all(abs(idx0 - selected) >= minSep)
        selected(end+1,1) = idx0; %#ok<SAGROW>
    end
end

% Sort selected minima by depth
[~, sidx] = sort(modC_smooth(selected), 'ascend');
selected = selected(sidx);

% Keep top few
nKeep = min(nKeep, length(selected));
cand_idx_ord = selected(1:nKeep);
cand_vals = modC_smooth(cand_idx_ord);

% Map back to time order
cand_idx_time = ord(cand_idx_ord);
cand_pts = Ybranch(cand_idx_time,:);

% ------------------------------------------------------------
% STEP 1: direct period-3 test on each candidate
% ------------------------------------------------------------
err3 = NaN(nKeep,1);

for i = 1:nKeep
    n = cand_idx_time(i);
    if n+3 <= size(Ybranch,1)
        err3(i) = norm(Ybranch(n+3,:) - Ybranch(n,:));
    end
end

% Candidate orbit points = small 3-step return error
orbit_ids = find(err3 < orbitTol);

% If too few found, take the best few automatically
if length(orbit_ids) < 3
    [~, idxsort] = sort(err3, 'ascend');
    orbit_ids = idxsort(1:min(3,length(idxsort)));
end

% ------------------------------------------------------------
% STEP 2: among orbit candidates, look for a genuine 3-cycle
% by nearest image matching within the orbit set
% ------------------------------------------------------------
orbit_cycle = [];
if length(orbit_ids) >= 3
    orbit_pts = cand_pts(orbit_ids,:);
    orbit_time = cand_idx_time(orbit_ids);

    next_map = NaN(length(orbit_ids),1);
    for a = 1:length(orbit_ids)
        n = orbit_time(a);
        if n+1 <= size(Ybranch,1)
            y_next = Ybranch(n+1,:);
            d = vecnorm(orbit_pts - y_next, 2, 2);
            [~, j] = min(d);
            next_map(a) = j;
        end
    end

    for a = 1:length(orbit_ids)
        b = next_map(a);
        if isnan(b), continue; end
        c = next_map(b);
        if isnan(c), continue; end
        d = next_map(c);
        if d == a && a~=b && b~=c && c~=a
            orbit_cycle = orbit_ids([a b c]);
            break
        end
    end
end

% If no exact symbolic 3-cycle found, use the 3 best return-error minima
if isempty(orbit_cycle)
    [~, idxsort] = sort(err3, 'ascend');
    orbit_cycle = idxsort(1:min(3,length(idxsort)));
end

% ------------------------------------------------------------
% STEP 3: preimage candidates
% small 1-step distance into orbit, but not small 3-step return
% ------------------------------------------------------------
orbit_pts = cand_pts(orbit_cycle,:);
preimage_ids = [];

for i = 1:nKeep
    if ismember(i, orbit_cycle)
        continue
    end

    n = cand_idx_time(i);
    if n+1 <= size(Ybranch,1)
        y_next = Ybranch(n+1,:);
        d1 = min(vecnorm(orbit_pts - y_next, 2, 2));

        if d1 < preTol && (isnan(err3(i)) || err3(i) >= orbitTol)
            preimage_ids(end+1,1) = i; %#ok<SAGROW>
        end
    end
end

% ------------------------------------------------------------
% REPORT
% ------------------------------------------------------------
fprintf('\nCandidate minima table:\n');
fprintf('candID   timeIndex   smoothedLogModulus   err3\n');
for i = 1:nKeep
    fprintf('%5d   %9d   %18.6f   %10.6g\n', ...
        i, cand_idx_time(i), cand_vals(i), err3(i));
end

fprintf('\nChosen orbit-cycle candidate IDs:\n');
disp(orbit_cycle(:)')

fprintf('Chosen orbit-cycle time indices:\n');
disp(cand_idx_time(orbit_cycle)')

fprintf('Candidate preimage IDs:\n');
disp(preimage_ids(:)')

if ~isempty(preimage_ids)
    fprintf('Candidate preimage time indices:\n');
    disp(cand_idx_time(preimage_ids)')
end

% ------------------------------------------------------------
% PLOTS
% ------------------------------------------------------------
figure(100)
plot(modC_ord, 'Color', [0.75 0.75 0.75]); hold on
plot(modC_smooth, 'b', 'LineWidth', 1.5)
scatter(cand_idx_ord, modC_smooth(cand_idx_ord), 50, 'k', 'filled')
scatter(cand_idx_ord(orbit_cycle), modC_smooth(cand_idx_ord(orbit_cycle)), 80, 'r', 'filled')
if ~isempty(preimage_ids)
    scatter(cand_idx_ord(preimage_ids), modC_smooth(cand_idx_ord(preimage_ids)), 80, 'b', 'filled')
end
xlabel('Ordered point index')
ylabel('log|\psi|')
title('Modulus minima: orbit candidates (red), preimages (blue)')
grid on

figure(101)
hold on
scatter(Ybranch(:,1), Ybranch(:,2), 8, [0.85 0.85 0.85], 'filled')
scatter(cand_pts(orbit_cycle,1), cand_pts(orbit_cycle,2), 90, 'r', 'filled')
if ~isempty(preimage_ids)
    scatter(cand_pts(preimage_ids,1), cand_pts(preimage_ids,2), 90, 'b', 'filled')
end
xlabel('PC1')
ylabel('PC2')
title('Candidate period-3 orbit (red) and preimages (blue)')
axis equal
grid on
legend('branch','period-3 candidates','preimage candidates','Location','best')








%% ============================================================
% PARTITION THE BRANCH USING THE DETECTED PERIOD-3 ORBIT
% AND BUILD A 3x3 TRANSITION MATRIX FOR THE SECOND-RETURN MAP
%
% REQUIRED VARIABLES ALREADY IN WORKSPACE:
%   Ybranch        : N x 2 points on one selected branch
%   ord            : ordering of branch points along the curve
%   orbit_cycle    : 1 x 3 indices of the detected period-3 orbit
%   cand_idx_ord   : ordered-index locations of candidate minima
%   cand_idx_time  : time-order locations of candidate minima
%
% OUTPUTS:
%   labels_ord_3   : region labels (1,2,3) in ordered-branch coordinates
%   labels_time_3  : same labels in time order
%   P3             : 3 x 3 transition matrix under second-return map
% ============================================================

% --- Get the ordered locations of the 3 orbit points ---
orbit_pos_ord = sort(cand_idx_ord(orbit_cycle));

% Sanity check
disp('Ordered branch locations of period-3 orbit points:')
disp(orbit_pos_ord(:)')

N = length(ord);

% --- Build 3 interval labels in ordered coordinates ---
% We cut the ordered branch at the three orbit locations.
% The intervals are:
%   Region 1: [orbit1, orbit2)
%   Region 2: [orbit2, orbit3)
%   Region 3: the wrap-around remainder
labels_ord_3 = zeros(N,1);

i1 = orbit_pos_ord(1);
i2 = orbit_pos_ord(2);
i3 = orbit_pos_ord(3);

labels_ord_3(i1:i2-1) = 1;
labels_ord_3(i2:i3-1) = 2;
labels_ord_3(i3:N)    = 3;
labels_ord_3(1:i1-1)  = 3;   % wrap-around part

% --- Map labels back to time order ---
labels_time_3 = zeros(N,1);
labels_time_3(ord) = labels_ord_3;

% --- Compute transition matrix for second-return map ---
P3 = zeros(3,3);

for n = 1:N-1
    a = labels_time_3(n);
    b = labels_time_3(n+1);
    P3(a,b) = P3(a,b) + 1;
end

% Normalize rows
rowSums = sum(P3,2);
for i = 1:3
    if rowSums(i) > 0
        P3(i,:) = P3(i,:) / rowSums(i);
    end
end

disp('3-state transition matrix P3:')
disp(P3)

% --- Plot the ordered branch with the 3 regions ---
figure
scatter(Ybranch(ord,1), Ybranch(ord,2), 16, labels_ord_3, 'filled')
xlabel('PC1')
ylabel('PC2')
title('3-region partition induced by period-3 orbit')
axis equal
grid on
box on
colormap(lines(3))
cb = colorbar;
cb.Ticks = [1 2 3];
cb.TickLabels = {'R1','R2','R3'};

% Highlight the orbit points
hold on
scatter(Ybranch(ord(orbit_pos_ord),1), Ybranch(ord(orbit_pos_ord),2), ...
    80, 'k', 'filled')

% --- Plot symbolic itinerary in time order ---
figure
stairs(labels_time_3(1:min(400,N)), 'LineWidth', 1.2)
xlabel('Iteration of second-return map')
ylabel('Region label')
title('3-region symbolic itinerary')
ylim([0.5 3.5])
yticks([1 2 3])
grid on

% --- Also show the labels along the ordered branch index ---
figure
stairs(labels_ord_3, 'LineWidth', 1.2)
xlabel('Ordered branch index')
ylabel('Region label')
title('3-region partition along ordered branch')
ylim([0.5 3.5])
yticks([1 2 3])
grid on

% --- Dominant next-state summary ---
[vals, idx] = max(P3, [], 2);
disp('Most likely next state from each current region:')
for i = 1:3
    fprintf('Region %d -> Region %d (probability %.3f)\n', i, idx(i), vals(i));
end