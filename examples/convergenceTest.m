function convergenceTest()
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
    mrstModule add vem mpfa mpsaw vemmech

    %% params setting
    % nref     : degree of refinement
    nref = 3; % default setting
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
                        'alpha'   , 1, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g1_k1_a1_eta0';
        savethisfigure(filename);
    end

    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 1, ...
                        'gridtype', 2, ... % Triangular grid, 90 degree angles
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g2_k1_a1_eta3';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 1, ...
                        'gridtype', 3, ... % Triangular grid, 90 degree angles
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g3_k1_a1_eta3';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 1, ...
                        'gridtype', 4, ... % Triangular grid, 90 degree angles
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g4_k1_a1_eta3';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 3, ...
                        'alpha'   , 1, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g1_k3_a1_eta3';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd2_g1_k1_a10_eta0';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 1, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd3_g1_k1_a1_eta0';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 3, ...
                        'alpha'   , 1, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd3_g1_k3_a1_eta0';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd3_g1_k1_a10_eta0';
        savethisfigure(filename);
    end
    
    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 1, ...
                        'gridtype', 5, ... % Cartesian Grid
                        'eta'     , 1e-8);
        
        output = mpsaPaperConvergenceFunc(params, 'doVem', doVem);
        figure
        plotConv(output, params);
        filename = 'convd3_g5_k1_a1_eta0';
        savethisfigure(filename);
    end
end

function plotConv(output, params)
    
    deL2 = output.deL2;
    
    Nd = params.Nd;
    nref = params.nref;
    kappa = params.kappa;
    alpha = params.alpha;
    gridtype = params.gridtype;
    eta = params.eta;

    log2N = (1 : nref)';
    plot(log2N, log2(deL2), '*-', 'linewidth', 4);
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel('log2(error)', 'fontsize', 18);
    xlabel('-log2(1/N)', 'fontsize', 18);
    caseTitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha);
    title(caseTitle);    
end


function casetitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha)
    dimstr = sprintf('%dD', Nd);
    switch gridtype
      case 1
        gridtypestr = 'Cartesian grid';
      case 2
        gridtypestr = 'Triangular grid, 90 degree angles';
      case 3
        gridtypestr = 'Equilateral triangles';
    end
    
    casetitle = sprintf('%s - %s, \\kappa = %0.3g, \\alpha = %0.3g, \\eta = %0.3g', dimstr, ...
                        gridtypestr, kappa, alpha, eta);
        
end

function savethisfigure(filename)
    savedir = '/home/xavier/Dropbox/figs/';
    filename = sprintf('%s.png', filename);
    filename = fullfile(savedir, filename);
    saveas(gcf, filename);
end
