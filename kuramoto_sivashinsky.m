% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE KURAMOTO-SIVASHINSKY POINCARE MAP
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script studies the Poincare map associated with a Galerkin
% projection of the Kuramoto-Sivashinsky equation at viscosity
% parameter nu = 0.0298.
%
% The workflow is:
%   1. Simulate the Galerkin system and extract Poincare section data
%   2. Load previously saved section data for repeated analysis
%   3. Compute a principal component representation of the section
%   4. Visualize the singular value spectrum and projected section
%   5. Build delay-coordinate evaluation matrices
%   6. Apply riggedDMD to approximate spectral measures and a
%      generalized eigenfunction
%   7. Visualize phase, ordering, and eigenfunction modulus
%
% In this script, the observable used for riggedDMD is the first
% principal component coordinate of the Poincare section data.
%
% MODEL
% The Galerkin system is a finite-dimensional truncation of the
% Kuramoto-Sivashinsky equation with m modes and viscosity nu.
%
% MAIN USER PARAMETERS
%   nu              : viscosity parameter
%   m               : number of Galerkin modes
%   dt              : timestep for ODE integration
%   dim             : number of delays in Hankel embedding
%   epsilon         : smoothing parameter for riggedDMD
%   order           : kernel order
%   thetaTarget     : angle for generalized eigenfunction
%   thetaGrid       : angle grid for spectral measure
%
% NOTES
%   - For repeated Koopman computations, the script loads saved
%     Poincare section data rather than re-simulating each time.
%   - The Poincare section data is projected onto principal components.
%   - If y1, y2, y3 denote the principal component coordinates, then one
%     may form the branch data via
%
%         Ybranch = [y1 y2 y3];
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

% Kuramoto-Sivashinsky parameters
nu = 0.0298;
m = 32;

% Time integration parameters
dt = 1e-3;
nSteps = 1e5;
tspan = (0:nSteps-1) * dt;

% ODE solver tolerances
odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Initial condition, chosen near a period-1 UPO
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
   -0.02766516
];
x0(15:m) = zeros(m-14,1);

% riggedDMD parameters
dim = 20;                    % number of delays
epsilon = 0.75;               % smoothing parameter
order = 2;                   % kernel order
thetaTarget = 0.25*pi;       % generalized eigenfunction angle
thetaGrid = -pi:0.01:pi;     % angles for spectral measure

%% Simulate the Galerkin system

[t, sol] = ode45(@(t,x) ks_rhs(x,nu,m), tspan, x0, odeOptions); %#ok<ASGLU>

%% Plot trajectory in the (a1,a2)-plane

figure(1)
clf
hold on

% Plot a portion of the attractor
plot(sol(1:50000,1), sol(1:50000,2), 'k', 'LineWidth', 2)

% Plot the Poincare section a2 = 0
a1Line = -1.5:0.01:1.5;
plot(a1Line, 0*a1Line, 'Color', [0 120/255 0], 'LineWidth', 5)

set(gca,'FontSize',16)
xlabel('$a_1$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$a_2$','Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_a1a2')

%% Construct Poincare section data

% We record points when a2 crosses from negative to nonnegative.
% The retained coordinates are [a1, a3, ..., am].
poincareData = [];
count = 1;

for n = 1:size(sol,1)-1
    if sol(n,2) < 0 && sol(n+1,2) >= 0
        poincareData(count,:) = sol(n+1,[1 3:m]); %#ok<SAGROW>
        count = count + 1;
    end
end

% Remove a transient crossing
poincareData = poincareData(2:end,:);

%% Load saved Poincare section data for repeated analysis

close all
clc

% This file is assumed to contain the variable x holding the section data.
load('ks_psec_0298_x2=0.mat')

% Rename for clarity
poincareData = x;

%% Compute principal component coordinates

% Use every second point, consistent with the second-return map analysis
dataForSVD = poincareData(2:2:end,:);

[U, S, V] = svd(dataForSVD, 'econ');

% Principal component coordinates
y1 = U(:,1) * S(1,1);
y2 = U(:,2) * S(2,2);
y3 = U(:,3) * S(3,3);

% Optional combined branch variable
Ybranch = [y1 y2 y3];

%% Plot singular value spectrum

figure(2)
clf

singularWeights = diag(S) / sum(diag(S));
bar(singularWeights, 'FaceColor', [36/255 122/255 254/255], 'EdgeColor', 'none')
hold on

% Mark very small singular values
smallThreshold = 1e-3;
isSmall = singularWeights < smallThreshold;
plot(find(isSmall), singularWeights(isSmall), 'k*', 'MarkerSize', 8, 'LineWidth', 1.5)

set(gca,'FontSize',16,'Box','on')
xlabel('Index','FontSize',16)
ylabel('Normalized singular values','FontSize',16)
grid on

save_figure_pdf(gcf, 'ks_singular_values')

%% Plot projected Poincare section and organizing orbits

addpath('ks_upos/')

figure(3)
clf
hold on

% Projected Poincare section data
plot(y1, y2, 'k.', 'MarkerSize', 10)

% Load period-3 orbit data
load('per3asym.mat')

% Project the two period-3 orbits into the first two PCs
per3proj1 = per3asymPsec1 * V(1:27,1:2);
per3proj2 = per3asymPsec2 * V(1:27,1:2);

% Plot period-3 orbits
plot(per3proj1(:,1), per3proj1(:,2), 'diamond', ...
    'MarkerSize', 16, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

plot(per3proj2(:,1), per3proj2(:,2), 'square', ...
    'MarkerSize', 16, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

% Plot a preimage point
per3Preimage = [8.15; -1.98]; %[8.2498; -1.9837];
plot(per3Preimage(1), per3Preimage(2), 'pentagram', ...
    'MarkerSize', 20, ...
    'Color', [151/255 1/255 200/255], ...
    'MarkerFaceColor', [151/255 14/255 200/255])

set(gca,'FontSize',16)
xlabel('PC$_1$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('PC$_2$','Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on

% Labels for the red orbit
text(per3proj2(1,1)-0.2,  per3proj2(1,2)+0.3, 'C', 'FontSize', 20)
text(per3proj2(2,1)-0.2,  per3proj2(2,2)+0.3, 'A', 'FontSize', 20)
text(per3proj2(3,1)-0.2,  per3proj2(3,2)+0.3, 'B', 'FontSize', 20)

% Labels for the green orbit
text(per3proj1(1,1)-0.15, per3proj1(1,2)+0.27, 'F', 'FontSize', 20)
text(per3proj1(2,1)-0.18, per3proj1(2,2)+0.27, 'D', 'FontSize', 20)
text(per3proj1(3,1)+0.05, per3proj1(3,2)-0.27, 'E', 'FontSize', 20)

% Label for the preimage
text(per3Preimage(1)+0.05, per3Preimage(2)-0.27, 'G', 'FontSize', 20)

save_figure_pdf(gcf, 'ks_psec')

%% Build Hankel evaluation matrices

observable = y1 - mean(y1);
observable = observable / std(observable);

H = hankel(observable);
PsiX = H(1:end-dim, 1:dim);
PsiY = H(1:end-dim, 2:dim+1);

%% Apply riggedDMD

addpath('main_routines/')

G = (PsiX' * PsiX) / size(PsiX,1);

coeffs = G / G(1,1);
coeffs = conj(coeffs(:));
coeffs = [coeffs; zeros(dim-2,1)];

G = toeplitz(coeffs, coeffs');
A = G(1:dim, 2:dim+1); 
G = G(1:dim, 1:dim);   

gCoeffs = zeros(dim,1);
gCoeffs(1) = 1;

[gModes, spectralDensity] = riggedDMD( ...
    PsiX, PsiY, 1/length(poincareData), epsilon, thetaTarget, [], ...
    'order', order, ...
    'g_coeffs', gCoeffs, ...
    'TH2', thetaGrid);

gModes = squeeze(gModes);

% Evaluate generalized eigenfunction
C = PsiX * gModes;

%% Plot spectral measure

figure(4)
clf
plot(thetaGrid, spectralDensity, 'Color', [1 69/255 79/255], 'LineWidth', 4)

xlabel('$\theta$','Interpreter','latex','FontSize',24,'FontWeight','bold')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
set(gca,'FontSize',16)
xlim([0 pi])
grid on

save_figure_pdf(gcf, 'ks_spec_measure')

%% Compute an ordering of the projected section and unwrap phase

theta = angle(C);

% Restrict to the points on which C is defined
P = [y1(1:length(theta)), y2(1:length(theta))];

% Order points along the branch
ord = greedy_curve_order(P);

thetaOrdered = unwrap(theta(ord));
modulusOrdered = log(abs(C(ord))); 

figure(5)
clf
plot(thetaOrdered, '.-', 'Color', [36/255 122/255 254/255], 'LineWidth', 2)
xlabel('Ordered data point index','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('Unwrapped phase','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
grid on
xlim([1 length(thetaOrdered)])

save_figure_pdf(gcf, 'ks_unwrapped')

%% Visualize the ordering along the attractor

figure(6)
clf
scatter(P(:,1), P(:,2), 12, 1:length(ord), 'filled')
colorbar
set(gca,'FontSize',16)
xlabel('PC$_1$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('PC$_2$','Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on

%% Plot generalized eigenfunction modulus

xPlot = y1(1:end-dim);
[~, sortIdx] = sort(xPlot);

figure(7)
clf
hold on

% Red period-3 orbit
xline(per3proj2(1,1), '--', 'Color', [1 69/255 79/255], 'LineWidth', 2)
xline(per3proj2(2,1), '--', 'Color', [1 69/255 79/255], 'LineWidth', 2)
xline(per3proj2(3,1), '--', 'Color', [1 69/255 79/255], 'LineWidth', 2)

% Green period-3 orbit
xline(per3proj1(1,1), '--', 'Color', [0 120/255 0], 'LineWidth', 2)
xline(per3proj1(2,1), '--', 'Color', [0 120/255 0], 'LineWidth', 2)
xline(per3proj1(3,1), '--', 'Color', [0 120/255 0], 'LineWidth', 2)

% Preimage point
xline(per3Preimage(1), '--', 'Color', [151/255 1/255 200/255], 'LineWidth', 2)

% Modulus plot
plot(xPlot(sortIdx), log(abs(C(sortIdx))), 'k-', 'LineWidth', 1.5)

xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\log|\varphi(x)|$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
axis([7.75 10.6 -5.5 0])
box on
grid on

save_figure_pdf(gcf, 'ks_eigenfunction')

%% Local functions

function dx = ks_rhs(x, nu, modes)
    %KS_RHS Right-hand side for a Galerkin truncation of KSE.
    %
    %   dx = ks_rhs(x, nu, modes)
    %
    % Inputs:
    %   x      - current state vector
    %   nu     - viscosity parameter
    %   modes  - number of Galerkin modes
    %
    % Output:
    %   dx     - time derivative
    
    dx = zeros(modes,1);
    
    for k = 1:modes
        dx(k) = (k^2) * (1 - nu*k^2) * x(k);
    
        for n = 1:(modes-k)
            dx(k) = dx(k) + 0.5 * k * x(n) * x(n+k);
        end
    
        for j = 1:(k-1)
            dx(k) = dx(k) - 0.25 * k * x(j) * x(k-j);
        end
    end
end

function ord = greedy_curve_order(P)
    %GREEDY_CURVE_ORDER Construct a nearest-neighbour ordering of point cloud.
    %
    %   ord = greedy_curve_order(P)
    %
    % Input:
    %   P    - N x 2 array of planar points
    %
    % Output:
    %   ord  - ordering of the points obtained by a greedy nearest-neighbour
    %          traversal starting from the leftmost point
    
    N = size(P,1);
    ord = zeros(N,1);
    used = false(N,1);
    
    % Start from the leftmost point
    [~, startIdx] = min(P(:,1));
    ord(1) = startIdx;
    used(startIdx) = true;
    
    for k = 2:N
        lastPoint = P(ord(k-1), :);
        dist2 = sum((P - lastPoint).^2, 2);
        dist2(used) = inf;
        [~, nextIdx] = min(dist2);
        ord(k) = nextIdx;
        used(nextIdx) = true;
    end
end

function save_figure_pdf(figHandle, fileName)
    %SAVE_FIGURE_PDF Save a figure to PDF using its on-screen size.
    
    set(figHandle, 'PaperUnits', 'centimeters');
    set(figHandle, 'Units', 'centimeters');
    
    pos = get(figHandle, 'Position');
    
    set(figHandle, 'PaperSize', [pos(3) pos(4)]);
    set(figHandle, 'PaperPositionMode', 'manual');
    set(figHandle, 'PaperPosition', [0 0 pos(3) pos(4)]);
    
    print(figHandle, '-dpdf', fileName);
end