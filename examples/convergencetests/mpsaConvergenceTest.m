mrstModule add vem mpfa mpsaw vemmech libgeometry

clear all
close all

%% params setting
% nref     : degree of refinement
% Nd       : dimension (2D or 3D)
% kappa    : Value of heterogenity in the domain (see paper)
% alpha    : Parameter to setup Lamé coefficient lambda (lambda = alpha*mu)
% gridtype : grid type (see mpsaPaperConvergenceFunc)
% eta      : Value used to set the position of the continuity point

% Possibility to run vem for comparison
doVem = false;

%% New Case
dothiscase = true;
params = struct('nref'    , 4, ...
                'Nd'      , 3, ...
                'kappa'   , 10, ...
                'alpha'   , 0, ...
                'gridtype', 5, ... % Cartesian
                'eta'     , 1/4);

output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', 100);

figure
hold on
plotConvTest(output, params);

