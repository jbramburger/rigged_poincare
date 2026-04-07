% -------------------------------------------------------------------------
% RIGGEDDMD ANALYSIS OF THE ROSSLER POINCARE MAP
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script simulates the Rossler system, constructs a Poincare section
% at x = 0 using crossings with x changing from negative to nonnegative,
% and applies rigged Dynamic Mode Decomposition (riggedDMD) to analyze the
% induced return dynamics.
%
% The script produces:
%   1. The Poincare section in the (y,z)-plane
%   2. The one-dimensional return map y_n -> y_{n+1}
%   3. A histogram approximating the invariant density in y
%   4. A riggedDMD approximation of the spectral measure
%   5. A generalized eigenfunction associated with a chosen angle
%
% MODEL
% The Rossler system is
%
%   x' = -y - z
%   y' =  x + a y
%   z' =  b + z(x-c)
%
% with parameter values
%   a = 0.1,  b = 0.1,  c = 18.
%
% DATA-DRIVEN SETUP
% The observable is taken to be the y-coordinate of the Poincare section,
% centered and normalized before constructing a Hankel delay embedding.
% The resulting dictionary evaluation matrices are then used in a riggedDMD
% computation.
%
% RIGGEDDMD IMPLEMENTATION
% The riggedDMD code used here is taken from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% MAIN USER PARAMETERS
%   dt              : ODE timestep
%   dim             : number of delays in Hankel embedding
%   epsilon         : smoothing parameter
%   order           : kernel order
%   thetaTarget     : angle for generalized eigenfunction
%   thetaGrid       : angle grid for spectral measure
%
% NOTES
%   - All Koopman computations are fully data-driven.
%   - The numerical integration uses a fixed timestep of 1e-3.
%   - For repeated riggedDMD experiments, it is often convenient to save
%     the Poincare section data and reload it rather than re-simulating.
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

% Rossler parameters
a = 0.1;
b = 0.1;
c = 18;

% Time integration parameters
dt = 1e-3;
nSteps = 1e8;
tspan = (0:nSteps-1) * dt;
x0 = [0.1; -8; 0.3];

% ODE solver tolerances
odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

% riggedDMD parameters
dim = 50;                       % number of delays in Hankel embedding
epsilon = 0.9;                  % smoothing parameter
order = 2;                      % kernel order
thetaTarget = 2*pi/3;           % angle for generalized eigenfunction
thetaGrid = -pi:0.01:pi;        % angles for spectral measure

% Known dynamically significant points in the return map
fixedPointY = -22.9130;
fixedPointPreimageY = -16.8105;

%% Simulate the Rossler system

[~, trajectory] = ode45(@(t,x) rossler_rhs(x,a,b,c), tspan, x0, odeOptions);

%% Construct Poincare section at x = 0 with positive crossing

% We record (y,z) whenever x crosses from negative to nonnegative.
poincareYZ = [];
count = 1;

for n = 1:size(trajectory,1)-1
    if trajectory(n,1) < 0 && trajectory(n+1,1) >= 0
        poincareYZ(count,:) = trajectory(n+1,2:3); %#ok<SAGROW>
        count = count + 1;
    end
end

% Remove an initial transient crossing
poincareYZ = poincareYZ(2:end,:);

%% Plot Poincare section in the (y,z)-plane

figure(1)
clf
plot(poincareYZ(:,1), poincareYZ(:,2), 'k.', 'MarkerSize', 10)
set(gca,'FontSize',16)
xlabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$z$','Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on
axis([-28 -11 0.00505 0.0055])

save_figure_pdf(gcf, 'rossler_2D_psec')

%% Plot one-dimensional return map in y

yData = poincareYZ(:,1);

figure(2)
clf
hold on

% Diagonal y_{n+1} = y_n
plot(yData(1:end-1), yData(1:end-1), '.', ...
    'Color', [0.5 0.5 0.5], 'MarkerSize', 5)

% Return map y_n -> y_{n+1}
plot(yData(1:end-1), yData(2:end), 'k.', 'MarkerSize', 10)

% Fixed point
plot(fixedPointY, fixedPointY, 'square', ...
    'MarkerSize', 20, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

% Preimage of fixed point
plot(fixedPointPreimageY, fixedPointY, 'diamond', ...
    'MarkerSize', 20, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

set(gca,'FontSize',16)
xlabel('$y_n$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$y_{n+1}$','Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on
box on
axis([-27 -11.5 -27 -11.5])

save_figure_pdf(gcf, 'rossler_1D_psec')

%% Plot histogram approximating the invariant density

figure(3)
clf
histogram(yData, 50, ...
    'Normalization', 'pdf', ...
    'FaceColor', [36/255 122/255 254/255])

set(gca,'FontSize',16)
xlabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('Density','Interpreter','latex','FontSize',24,'FontWeight','bold')
box on

save_figure_pdf(gcf, 'rossler_density')

%% Prepare data for riggedDMD

% For repeated experiments, one may prefer to load previously saved
% Poincare data rather than regenerate it every time.
close all
clc

% This file is assumed to contain the y-component of the Poincare data.
load rossler_psec.mat

% Center and normalize the observable
g = x - mean(x);
g = g / std(g);

% Build Hankel matrix and dictionary evaluation matrices
H = hankel(g);
PsiX = H(1:end-dim, 1:dim);
PsiY = H(1:end-dim, 2:dim+1);

%% Apply riggedDMD

addpath('main_routines/')

% Gram matrix
G = (PsiX' * PsiX) / size(PsiX,1);

% Toeplitz symmetrization used in the riggedDMD setup
coeffs = G / G(1,1);
coeffs = conj(coeffs(:));
coeffs = [coeffs; zeros(dim-2,1)];

G = toeplitz(coeffs, coeffs');
A = G(1:dim, 2:dim+1); 
G = G(1:dim, 1:dim);   

% Observable coefficients in the dictionary basis
gCoeffs = zeros(dim,1);
gCoeffs(1) = 1;

% Weighted riggedDMD
[gModes, spectralDensity] = riggedDMD( ...
    PsiX, PsiY, 1/length(x), epsilon, thetaTarget, [], ...
    'order', order, ...
    'g_coeffs', gCoeffs, ...
    'TH2', thetaGrid);

gModes = squeeze(gModes);

%% Plot spectral measure

figure(4)
clf
plot(thetaGrid, spectralDensity, ...
    'Color', [1 69/255 79/255], 'LineWidth', 4)

xlabel('$\theta$','Interpreter','latex','FontSize',24,'FontWeight','bold')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
set(gca,'FontSize',16)
xlim([0 pi])
grid on

save_figure_pdf(gcf, 'rossler_spec_measure')

%% Plot generalized eigenfunction modulus

% Evaluate generalized eigenfunction on the data
C = PsiX * gModes;

% Sort by the one-dimensional observable for visualization
xPlot = x(1:end-dim);
[~, sortIdx] = sort(xPlot);

figure(5)
clf
hold on

% Fixed point and preimage markers
xline(fixedPointY, '--', ...
    'Color', [1 69/255 79/255], 'LineWidth', 2)
xline(fixedPointPreimageY, '--', ...
    'Color', [0 120/255 0], 'LineWidth', 2)

% Modulus of generalized eigenfunction
plot(xPlot(sortIdx), log(abs(C(sortIdx))), 'k-', 'LineWidth', 1.5)

xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$\log|\varphi(x)|$','Interpreter','latex','FontSize',24,'FontWeight','bold')
set(gca,'FontSize',16)
axis([-27 -11.5 -7 0])
box on
grid on

save_figure_pdf(gcf, 'rossler_eigenfunction')

%% Local functions

function dx = rossler_rhs(x,a,b,c)
    %ROSSLER_RHS Right-hand side of the Rossler system.
    %
    %   dx = rossler_rhs(x,a,b,c)
    %
    % Inputs:
    %   x  - current state vector [x; y; z]
    %   a, b, c - Rossler parameters
    %
    % Output:
    %   dx - time derivative
    
    dx = [
        -x(2) - x(3);
         x(1) + a*x(2);
         b + x(3)*(x(1) - c)
    ];
end

function save_figure_pdf(figHandle, fileName)
    %SAVE_FIGURE_PDF Save a figure to PDF using its on-screen size.
    %
    %   save_figure_pdf(figHandle, fileName)
    
    set(figHandle, 'PaperUnits', 'centimeters');
    set(figHandle, 'Units', 'centimeters');
    
    pos = get(figHandle, 'Position');
    
    set(figHandle, 'PaperSize', [pos(3) pos(4)]);
    set(figHandle, 'PaperPositionMode', 'manual');
    set(figHandle, 'PaperPosition', [0 0 pos(3) pos(4)]);
    
    print(figHandle, '-dpdf', fileName);
end