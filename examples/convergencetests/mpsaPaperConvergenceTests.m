function mpsaPaperConvergenceTests()
% Script for running convergence test
% see function mpsaPaperConvergenceFunc

% reference paper:
% 
% title={Finite volume methods for elasticity with weak symmetry},
% author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
% journal={International Journal for Numerical Methods in Engineering},
% volume={112},
% number={8},
% pages={939--962},
% year={2017},
% publisher={Wiley Online Library}

    %% Load necessary modules
    mrstModule add vem mpfa mpsaw vemmech libgeometry

    close all
    
    %% params setting
    % nref     : degree of refinement
    nref = 3;  % default setting
    % Nd       : dimension (2D or 3D)
    % kappa    : Value of heterogenity in the domain (see paper)
    % alpha    : Parameter to setup Lamé coefficient lambda (lambda = alpha*mu)
    % gridtype : grid type (see mpsaPaperConvergenceFunc)
    % eta      : Value used to set the position of the continuity point
    
    % Possibility to run vem for comparison
    doVem = false;

    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        hold on
        plotConv(output, params);
        
        savethisfigure(params);
    end

    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 2, ... % Triangular grid, 90 degree angles
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end

    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 3, ... % Triangular grid, equi - alternate
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end

    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 4, ... % Triangular grid equi - non alternate
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          500);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          500);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          500);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ... % tetras
                        'gridtype', 5, ... 
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          100);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ... % tetras
                        'gridtype', 5, ... 
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          100);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end    
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ... % tetras
                        'gridtype', 5, ... 
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', ...
                                          100);
        figure
        plotConv(output, params);
        
        savethisfigure(params);
    end    
        
end

function savethisfigure(params)
    
    d = params.Nd;
    kappa = params.kappa;
    alpha = params.alpha;
    gridtype = params.gridtype;
    eta = params.eta;
    if eta == 1/3
        eta = 3;
    end
    filename = sprintf('convd_d%d_g%d_k%1.g_a%d_eta%d.png', d, gridtype, kappa, alpha, eta);
    savedir = '/home/xavier/Dropbox/figs/';
    filename = fullfile(savedir, filename);
    saveas(gcf, filename);
end
