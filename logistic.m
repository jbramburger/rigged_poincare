% -------------------------------------------------------------------------
% EDMD, MPEDMD, AND RIGGEDDMD FOR THE LOGISTIC MAP
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script generates trajectory data from the logistic map and applies
% several data-driven Koopman operator methods to analyze its spectral and
% topological structure.
%
% The script computes:
%   1. Standard EDMD approximations
%   2. Measure-preserving EDMD (mpEDMD)
%   3. riggedDMD approximations of spectral measures and generalized
%      eigenfunctions
%
% The observable is constructed using a delay-coordinate embedding
% (Hankel matrix), giving a Krylov-type approximation of the Koopman
% operator.
%
% The riggedDMD implementation used here is taken from:
%   https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition
%
% MAIN USER PARAMETERS
%   r               : logistic map parameter
%   nTransient      : number of transient iterates discarded
%   nData           : number of retained iterates
%   dim             : number of delays in the Hankel embedding
%   epsilon         : smoothing parameter for riggedDMD
%   order           : order of the rational kernel
%   thetaTarget     : angle for generalized eigenfunction
%   thetaGrid       : angle grid for spectral measure
%
% NOTES
%   - The observable is centered and normalized before analysis.
%   - All computations are fully data-driven and use only trajectory data.
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear
close all
clc

%% User parameters

% Logistic map parameter
r = 3.7;

% Data generation parameters
nTransient = 100;
nData = 20000;
x0 = 0.25;

% Delay embedding parameter
dim = 51;

% riggedDMD parameters
thetaTarget = pi;
thetaGrid = -pi:0.01:pi;
epsilon = 0.01;
order = 2;

% mpEDMD eigenfunction index to visualize
eigenfunctionIndex = 23;

% Kneading sequence length
kneadingLength = 40;
criticalPoint = 0.5;

%% Generate trajectory data

nTotal = nTransient + nData;

x = zeros(nTotal + 1, 1);
x(1) = x0;

for n = 1:nTotal
    x(n+1) = logistic_map(x(n), r);
end

% Remove transients
x = x(nTransient+1:end);

%% Compute the kneading sequence

criticalOrbit = zeros(kneadingLength + 1, 1);
criticalOrbit(1) = criticalPoint;

for n = 1:kneadingLength
    criticalOrbit(n+1) = logistic_map(criticalOrbit(n), r);
end

kneadingSequence = char(zeros(kneadingLength,1));

for n = 2:kneadingLength+1
    if abs(criticalOrbit(n) - criticalPoint) < 1e-12
        kneadingSequence(n-1) = 'C';
    elseif criticalOrbit(n) < criticalPoint
        kneadingSequence(n-1) = 'L';
    else
        kneadingSequence(n-1) = 'R';
    end
end

disp('Kneading sequence:')
disp(kneadingSequence.')

%% Build dictionary evaluation matrices

observable = x - mean(x);
observable = observable / std(observable);

H = hankel(observable);
PsiX = H(1:end-dim, 1:dim);
PsiY = H(1:end-dim, 2:dim+1);

%% Apply standard EDMD

K_edmd = PsiX \ PsiY;
[V_edmd, D_edmd] = eig(K_edmd);

lambda_edmd = diag(D_edmd);
phi_edmd = PsiX * V_edmd; 

%% Apply mpEDMD

% QR-based construction
[Q, R] = qr(PsiX, "econ");
T = (R') \ (PsiY' * Q);

[U, ~, V] = svd(T);
[W_mp, D_mp] = schur(V * U', 'complex');

V_mp = R \ W_mp;
phi_mp = PsiX * V_mp;
lambda_mp = diag(D_mp);
K_mp = (R \ (V * U')) * R; 

%% Plot EDMD and mpEDMD eigenvalues

figure(1)
clf
hold on

% Unit circle
thetaCircle = linspace(0, 2*pi, 201);
xUnit = cos(thetaCircle);
yUnit = sin(thetaCircle);
plot(xUnit, yUnit, 'k--', 'LineWidth', 2)

% Coordinate axes
plot(1.1*xUnit, 0*yUnit, 'k', 'LineWidth', 1)
plot(0*xUnit, 1.1*yUnit, 'k', 'LineWidth', 1)

% EDMD eigenvalues
plot(real(lambda_edmd), imag(lambda_edmd), 'square', ...
    'MarkerSize', 16, ...
    'Color', [1 69/255 79/255], ...
    'MarkerFaceColor', [1 69/255 79/255])

% mpEDMD eigenvalues
plot(real(lambda_mp), imag(lambda_mp), 'diamond', ...
    'MarkerSize', 16, ...
    'Color', [0 120/255 0], ...
    'MarkerFaceColor', [0 120/255 0])

xlabel('$\mathrm{Re}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\mathrm{Im}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
axis([-1.1 1.1 -1.1 1.1])
grid on

save_figure_pdf(gcf, 'logistic_eigenvalues_2')

%% Plot an mpEDMD eigenfunction

xPlot = x(1:end-dim);
[~, sortIdx] = sort(xPlot);

figure(2)
clf
plot(xPlot(sortIdx), log(abs(phi_mp(sortIdx, eigenfunctionIndex))), ...
    'k-', 'LineWidth', 1.5)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\log|\varphi(x)|$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
axis([0.26 0.93 -14.5 -3.5])
grid on

save_figure_pdf(gcf, 'logistic_eigenfunction')

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
    PsiX, PsiY, 1/length(x), thetaOrDefault(epsilon), thetaTarget, [], ...
    'order', order, ...
    'g_coeffs', gCoeffs, ...
    'TH2', thetaGrid);

gModes = squeeze(gModes);

%% Plot spectral measure

figure(3)
clf
plot(thetaGrid, spectralDensity, ...
    'Color', [1 69/255 79/255], 'LineWidth', 4)

xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
set(gca, 'FontSize', 16)
xlim([0 pi])
grid on

save_figure_pdf(gcf, 'logistic_measure_2')

%% Plot generalized eigenfunction modulus

C = PsiX * gModes;

figure(4)
clf
plot(xPlot(sortIdx), log(abs(C(sortIdx))), 'k-', 'LineWidth', 1.5)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$\log|\varphi(x)|$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'FontSize', 16)
axis([0.26 0.93 -11 0.5])
grid on

save_figure_pdf(gcf, 'logistic_rigged_eigenfunction')

%% Local functions

function y = logistic_map(x, r)
    %LOGISTIC_MAP Evaluate the logistic map.
    %
    %   y = logistic_map(x, r)
    %
    % Inputs:
    %   x  - current state
    %   r  - logistic map parameter
    %
    % Output:
    %   y  - next iterate
    
    y = r * x .* (1 - x);
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

function epsOut = thetaOrDefault(epsIn)
    %THETAORDEFAULT Pass through the riggedDMD smoothing parameter.
    %
    % This helper is included only to keep the call structure visually
    % parallel with the other scripts.
    
    epsOut = epsIn;
end