function biotConvergenceTests()
% Script for running convergence test for MPSA-MPFA Biot
% see function biotConvergenceFunc

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

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    % Load necessary modules
    mrstModule add vem mpfa mpsaw vemmech libgeometry

    close all
    dosave = true;
    savecount = 1; % counter for setting up filenames.
    filerootname = 'biotconvoutput';
    
    % params is a struct with the fields:
    %    Nd       : Spatial dimension (Nd = 2 or 3)
    %    nref     : number of refinement level
    %    eta      : Value used to set the position of the continuity point
    %    gridtype : Different grid type can be used, (see function gridForConvTest)
    %               1: Cartesian
    %               2: Triangles by alternating bisection of triangles
    %               3: Equilateral triangles
    %               4: Triangles by uniform bisection
    %               5: Tetrehedral grid  
    %    mu       : Lame first coefficient
    %    lambda   : Lame second coefficcient
    %    K        : isotropic permeability, value of unique diagonal coefficient
    %    alpha    : Biot coefficient
    %    rho      : fluid weak compressibility coefficient
    %    tau      : time discretization coefficient (should be always set to one)

    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 2, ...
                        'gridtype', 1, ...
                        'eta'     , 0, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 2, ...
                        'gridtype', 2, ...
                        'eta'     , 1/3, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end        
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 2, ...
                        'gridtype', 3, ...
                        'eta'     , 1/3, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 2, ...
                        'gridtype', 4, ...
                        'eta'     , 1/3, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'gridtype', 1, ...
                        'eta'     , 0, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params, 'blocksize', 100);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'gridtype', 5, ...
                        'eta'     , 1/4, ...
                        'mu'      , 1, ...
                        'lambda'  , 1, ...
                        'alpha'   , 1, ...
                        'K'       , 1, ...
                        'tau'     , 1, ...
                        'rho'     , 1);
        
        output = biotConvergenceFunc(params, 'blocksize', 100);
        plotConv(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
    end
                                                            
end

function plotConv(output, params, varargin)
    
    opt = struct('savefigure', false);
    opt = merge_options(opt, varargin{:});
    
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
    caseTitle = setCaseTitle('displacement', params);
    axis equal
    title(caseTitle);    

    if opt.savefigure
        savethisfigure(params);
    end
    
    figure
    plot(log2N, log2(perrL2), '*-', 'linewidth', 4);
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel('log2(error)', 'fontsize', 18);
    xlabel('-log2(1/N)', 'fontsize', 18);
    title('error pressure');
    caseTitle = setCaseTitle('pressure', params);
    axis equal
    title(caseTitle);    
    
    if opt.savefigure
        savethisfigure(params);
    end
    
end


function casetitle = setCaseTitle(valstr, params)
    
    Nd = params.Nd;
    gridtype = params.gridtype;
    eta = params.eta;
    
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
    
    casetitle = sprintf('%s, %s - %s, \\eta = %0.3g', valstr, dimstr, gridtypestr, eta);
        
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
