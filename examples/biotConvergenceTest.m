function biotConvergenceTest()
% Script for running convergence test
% see function mpsaPaperConvergenceFunc

% reference papers:
% 
% title={Finite volume methods for elasticity with weak symmetry},
% author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
% journal={International Journal for Numerical Methods in Engineering},
% volume={112},
% number={8},
% pages={939--962},
% year={2017},
% publisher={Wiley Online Library}
%
% and
%
% title={Stable cell-centered finite volume discretization for Biot equations},
% author={Nordbotten, Jan Martin},
% journal={SIAM Journal on Numerical Analysis},
% volume={54},
% number={2},
% pages={942--968},
% year={2016},
% publisher={SIAM}


    %% Load necessary modules
    mrstModule add vem mpfa mpsaw vemmech libgeometry

    close all
    
    % params is a struct with the fields:
    %    Nd       : Spatial dimension (Nd = 2 or 3)
    %    nref     : number of refinement level
    %    eta      : Value used to set the position of the continuity point
    %    gridtype : Different grid type can be used, see below
    %               1 : Cartesian grid
    %               2 : Triangular grid, 90 degree angles
    %               3 : Equilateral triangles
    %    mu       : Lame first coefficient
    %    lambda   : Lame second coefficcient
    %    K        : isotropic permeability, value of unique diagonal coefficient
    %    alpha    : Biot coefficient
    %    rho      : fluid weak compressibility coefficient
    %    tau      : time discretization coefficient (should be always set to one)

    %% New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 2, ...
                        'gridtype', 1, ...
                        'eta'     , 0, ...
                        'mu'      , 1, ...
                        'lambda'  , 0, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 0);
        
        output = biotConvergenceFunc(params);
        plotConv(output, params);
        
        % savethisfigure(params);
    end
    

    %% New Case
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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
    dothiscase = false;
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

function plotConv(output, params)
    
    uerrL2 = output.uerrL2;
    perrL2 = output.perrL2;
    
    nref = params.nref;
    % Nd = params.Nd;
    % gridtype = params.gridtype;
    % eta = params.eta;

    figure
    log2N = (1 : nref)';
    plot(log2N, log2(uerrL2), '*-', 'linewidth', 4);
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel('log2(error)', 'fontsize', 18);
    xlabel('-log2(1/N)', 'fontsize', 18);
    title('error displacement');
    % caseTitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha);
    % title(caseTitle);    

    figure
    plot(log2N, log2(perrL2), '*-', 'linewidth', 4);
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel('log2(error)', 'fontsize', 18);
    xlabel('-log2(1/N)', 'fontsize', 18);
    title('error pressure');
    % caseTitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha);
    % title(caseTitle);    
    
end


function casetitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha)
    dimstr = sprintf('%dD', Nd);
    switch gridtype
      case 1
        gridtypestr = 'Cartesian grid';
      case 2
        gridtypestr = 'Triangular grid, 90 degree angles';
      case 3
        gridtypestr = 'Equilateral triangles (alt)';
      case 4
        gridtypestr = 'Equilateral triangles';
      case 5
        gridtypestr = 'Tetrahedrals';
    end
    
    casetitle = sprintf('%s - %s, \\kappa = %0.3g, \\alpha = %0.3g, \\eta = %0.3g', dimstr, ...
                        gridtypestr, kappa, alpha, eta);
        
end

function savethisfigure(params)
    error('do not want to save  now')
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
