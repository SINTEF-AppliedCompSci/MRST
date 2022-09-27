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
    filerootname = 'mpsaconvoutput';
    
    % params setting
    % nref     : degree of refinement
    % Nd       : dimension (2D or 3D)
    % kappa    : Value of heterogenity in the domain (see paper)
    % alpha    : Parameter to setup Lamé coefficient lambda (lambda = alpha*mu)
    % gridtype : grid type (see gridForConvTest)
    % eta      : Value used to set the position of the continuity point
    

    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        hold on
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end

    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 2, ... % Triangular grid, 90 degree angles
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end

    return

    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 3, ... % Triangular grid, equi - alternate
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end

    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 4, ... % Triangular grid equi - non alternate
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 6, ...
                        'Nd'      , 2, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 4, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', ...
                                          500);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', ...
                                          500);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ...
                        'gridtype', 1, ... % Cartesian Grid
                        'eta'     , 0);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', ...
                                          500);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 0, ... % tetras
                        'gridtype', 5, ... 
                        'eta'     , 1/3);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', ...
                                          100);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 10, ...
                        'alpha'   , 0, ... 
                        'gridtype', 5, ... % tetras
                        'eta'     , 1/4);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', 100);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
        savethisfigure(params);
    end    
    
    % New Case
    dothiscase = true;
    if dothiscase
        params = struct('nref'    , 5, ...
                        'Nd'      , 3, ...
                        'kappa'   , 1, ...
                        'alpha'   , 10, ... % tetras
                        'gridtype', 5, ... 
                        'eta'     , 1/4);
        
        output = mpsaPaperConvergenceFunc(params, 'blocksize', 100);
        figure
        plotConvTest(output, params);
        if dosave
            filename = sprintf('%s%d.mat', filerootname, savecount);
            save(filename, 'params', 'output');
            savecount = savecount + 1;
        end
        
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
    savedir = mrstOutputDirectory();
    filename = fullfile(savedir, filename);
    saveas(gcf, filename);
end
