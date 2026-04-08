% -------------------------------------------------------------------------
% PLOTTER FOR UNSTABLE PERIODIC ORBITS IN THE KURAMOTO-SIVASHINSKY EQUATION
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script visualizes selected unstable periodic orbits (UPOs) of a
% Galerkin truncation of the Kuramoto-Sivashinsky equation at viscosity
% parameter nu = 0.0298.
%
% The script produces plots in two complementary coordinate systems:
%
%   1. Physical Galerkin coordinates:
%      - UPOs are plotted in the (a1,a2)-plane
%      - A chaotic trajectory is shown in the background for context
%
%   2. Poincare-section coordinates:
%      - Poincare section data is projected onto the first two principal
%        components
%      - Period-1 intersections are plotted on top of the projected
%        section
%      - Period-3 organizing orbits and a distinguished preimage point
%        are included for reference
%
% INPUT FILES
% The script assumes the following files are available:
%   ks_psec_0298_x2=0.mat   - Poincare section data
%   per1sym.mat             - symmetric period-1 orbit and section point
%   per1asym.mat            - asymmetric period-1 orbit and section points
%   per3asym.mat            - asymmetric period-3 orbit section points
%
% NOTES
%   - The Poincare section data is projected using an SVD/PCA basis
%     computed from every second point, consistent with second-return
%     map analysis.
%   - If y1, y2, y3 are the principal component coordinates, then
%     one may form the branch data as
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
nModes = 32;

% Time integration parameters for background trajectory
dt = 1e-3;
nSteps = 1e5;
tspan = (0:nSteps-1) * dt;

% ODE solver tolerances
odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,nModes));

% Initial condition chosen near a period-1 UPO
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
x0(15:nModes) = zeros(nModes - 14, 1);

% Add directory containing UPO data
addpath('ks_upos/')

%% Load Poincare section data and compute principal components

load('ks_psec_0298_x2=0.mat')

% Rename for clarity
poincareData = x;

% Use every second point, consistent with second-return map analysis
dataForSVD = poincareData(2:2:end,:);

[U, S, V] = svd(dataForSVD, 'econ');

% Principal component coordinates
y1 = U(:,1) * S(1,1);
y2 = U(:,2) * S(2,2);
y3 = U(:,3) * S(3,3); 

% Optional combined branch variable
Ybranch = [y1 y2 y3]; 

%% Generate a chaotic trajectory for background plotting

[~, trajectory] = ode45(@(t,x) ks_rhs(x,nu,nModes), tspan, x0, odeOptions);

%% Load UPO data

% Period-1 data
load('per1sym.mat')
load('per1asym.mat')

% Period-2 data
load('per2asym.mat')

% Period-3 data
load('per3asym.mat')

%% Plot symmetric period-1 orbit in the (a1,a2)-plane

figure(1)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Symmetric period-1 orbit
plot([per1sym(:,1); per1sym(:,1)], [per1sym(:,2); per1sym(:,2)], ...
    'Color', [36/255 122/255 254/255], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per1_sym')

%% Plot asymmetric period-1 orbit in the (a1,a2)-plane

figure(2)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Asymmetric period-1 orbit
plot([per1asym1(:,1); per1asym1(:,1)], [per1asym1(:,2); per1asym1(:,2)], ...
    'Color', [0.85 0.33 0.0], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per1_asym')

%% Project UPO intersections into the principal component coordinates

% Period-1 intersections
per1projSym   = per1symPsec   * V(1:31,1:2);
per1projAsym1 = per1asymPsec1 * V(1:31,1:2);
per1projAsym2 = per1asymPsec2 * V(1:31,1:2);

% Period-2 intersections
per2proj1 = per2asymPsec1 * V(1:31,1:2); % 57 orbit
per2proj2 = per2asymPsec2 * V(1:31,1:2); % 23 orbit

% Period-3 intersections
% These can only be projected using the first 27 rows of V.
per3proj1 = per3asymPsec1 * V(1:27,1:2);
per3proj2 = per3asymPsec2 * V(1:27,1:2);
per3proj3 = per3asymPsec3 * V(1:27,1:2);
per3proj4 = per3asymPsec4 * V(1:27,1:2);

% Distinguished preimage point
per3Preimage = [8.15; -1.98];

%% Plot projected Poincare section with period-1 and period-3 data

figure(3)
clf
hold on

% Projected Poincare section data
plot(y1, y2, 'k.', 'MarkerSize', 10)

% Period-3 intersections
plot(per3proj1(:,1), per3proj1(:,2), 'diamond', ...
    'MarkerSize', 12, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

plot(per3proj2(:,1), per3proj2(:,2), 'square', ...
    'MarkerSize', 12, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

% Distinguished preimage
plot(per3Preimage(1), per3Preimage(2), 'pentagram', ...
    'MarkerSize', 14, ...
    'Color', [151/255 1/255 200/255], ...
    'MarkerFaceColor', [151/255 14/255 200/255])

% Period-1 intersections
plot(per1projSym(1), per1projSym(2), '.', ...
    'MarkerSize', 60, ...
    'Color', [36/255 122/255 254/255])

plot(per1projAsym1(1), per1projAsym1(2), '.', ...
    'MarkerSize', 60, ...
    'Color', [0.85 0.33 0.0])

plot(per1projAsym2(1), per1projAsym2(2), '.', ...
    'MarkerSize', 60, ...
    'Color', [0.85 0.33 0.0])

set(gca, 'FontSize', 16)
xlabel('PC$_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('PC$_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on

save_figure_pdf(gcf, 'ks_per1_psec')

%% Plot period-2 orbit in the (a1,a2)-plane

figure(4)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Period-2 orbit
plot([per2asym2(:,1); per2asym2(:,1)], [per2asym2(:,2); per2asym2(:,2)], ...
    'Color', [36/255 122/255 254/255], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per2_1')

%% Plot other period-2 orbit in the (a1,a2)-plane

figure(5)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Other period-2 orbit
plot([per2asym1(:,1); per2asym1(:,1)], [per2asym1(:,2); per2asym1(:,2)], ...
    'Color', [0.85 0.33 0.0], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per2_2')

%% Plot projected Poincare section with period-2 and period-3 data

figure(6)
clf
hold on

% Projected Poincare section data
plot(y1, y2, 'k.', 'MarkerSize', 10)

% Period-3 intersections
plot(per3proj1(:,1), per3proj1(:,2), 'diamond', ...
    'MarkerSize', 12, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

plot(per3proj2(:,1), per3proj2(:,2), 'square', ...
    'MarkerSize', 12, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

% Distinguished preimage
plot(per3Preimage(1), per3Preimage(2), 'pentagram', ...
    'MarkerSize', 14, ...
    'Color', [151/255 1/255 200/255], ...
    'MarkerFaceColor', [151/255 14/255 200/255])

% Period-2 intersections
plot(per2proj2(:,1), per2proj2(:,2), '.', ...
    'MarkerSize', 60, ...
    'Color', [36/255 122/255 254/255])

plot(per2proj1(:,1), per2proj1(:,2), '.', ...
    'MarkerSize', 60, ...
    'Color', [0.85 0.33 0.0])

set(gca, 'FontSize', 16)
xlabel('PC$_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('PC$_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on

save_figure_pdf(gcf, 'ks_per2_psec')

%% Plot first selected period-3 orbit in the (a1,a2)-plane

figure(7)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Period-3 orbit: per3asym1
plot([per3asym1(:,1); per3asym1(:,1)], [per3asym1(:,2); per3asym1(:,2)], ...
    'Color', [0 120/255 0], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per3_1')

%% Plot second selected period-3 orbit in the (a1,a2)-plane

figure(8)
clf
hold on

% Chaotic background trajectory
plot(trajectory(1:50000,1), trajectory(1:50000,2), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 2)

% Period-3 orbit: per3asym3
plot([per3asym3(:,1); per3asym3(:,1)], [per3asym3(:,2); per3asym3(:,2)], ...
    'Color', [0.85 0.33 0.0], 'LineWidth', 4)

set(gca, 'FontSize', 16)
xlabel('$a_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$a_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on
axis([-4 4 -1.5 1.5])

save_figure_pdf(gcf, 'ks_per3_2')

%% Plot projected Poincare section with selected period-3 intersections

figure(9)
clf
hold on

% Projected Poincare section data
plot(y1, y2, 'k.', 'MarkerSize', 10)

% Distinguished preimage
plot(per3Preimage(1), per3Preimage(2), 'pentagram', ...
    'MarkerSize', 14, ...
    'Color', [151/255 1/255 200/255], ...
    'MarkerFaceColor', [151/255 14/255 200/255])

% Period-3 intersections
plot(per3proj1(:,1), per3proj1(:,2), 'diamond', ...
    'MarkerSize', 16, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

plot(per3proj2(:,1), per3proj2(:,2), 'square', ...
    'MarkerSize', 16, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

plot(per3proj3(:,1), per3proj3(:,2), '.', ...
    'MarkerSize', 50, ...
    'Color', [0.85 0.33 0.0], ...
    'MarkerFaceColor', [0.85 0.33 0.0])

plot(per3proj4(:,1), per3proj4(:,2), '.', ...
    'MarkerSize', 50, ...
    'Color', [36/255 122/255 254/255], ...
    'MarkerFaceColor', [36/255 122/255 254/255])


set(gca, 'FontSize', 16)
xlabel('PC$_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('PC$_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
grid on
box on

save_figure_pdf(gcf, 'ks_per3_psec')



%% Local functions

function dx = ks_rhs(x, nu, nModes)
    %KS_RHS Right-hand side for a Galerkin truncation of KSE.
    %
    %   dx = ks_rhs(x, nu, nModes)
    %
    % Inputs:
    %   x       - current state vector
    %   nu      - viscosity parameter
    %   nModes  - number of Galerkin modes
    %
    % Output:
    %   dx      - time derivative
    
    dx = zeros(nModes,1);
    
    for k = 1:nModes
        dx(k) = (k^2) * (1 - nu*k^2) * x(k);
    
        for n = 1:(nModes-k)
            dx(k) = dx(k) + 0.5 * k * x(n) * x(n+k);
        end
    
        for j = 1:(k-1)
            dx(k) = dx(k) - 0.25 * k * x(j) * x(k-j);
        end
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