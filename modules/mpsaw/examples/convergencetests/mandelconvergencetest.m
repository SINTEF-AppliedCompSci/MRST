function outputs = mandelconvergencetest
%Undocumented Utility Function

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


    % Spatial convergence
    params = setDefaultParameters();
    params.rampup = 1;
    ts = params.fixedtsteps;
    nxs = [5; 10; 100];
    nnxs = numel(nxs);
    states = cell(nnxs, 1);
    
    close all
    
    outputs = {};
    
    for i = 1 : nnxs
        params.nx = nxs(i);
        output = mandelrun(params);
        outputs{end + 1} = output;
        figure
        hold on
        plotsol(output);
        plotexactsol(output, 1e-2);
        titlestr = sprintf('nx = %d, dt = %g', params.nx, params.dt);
        title(titlestr)
        legstr = arrayfun(@(x) sprintf('%g', x), ts, 'uniformoutput', false);
        h = legend(legstr, 'location', 'eastoutside');
    end
    
    
    % Time convergence
    params = setDefaultParameters();
    params.rampup = 1;
    dts = [1e-2; 1e-3];
    ndts = numel(dts);
    states = cell(ndts, 1);
    
    
    for i = 1 : ndts
        params.dt = dts(i);
        output = mandelrun(params);
        outputs{end + 1} = output;
        figure
        hold on
        plotsol(output);
        plotexactsol(output, 1e-2);
        titlestr = sprintf('nx = %d, dt = %g', params.nx, params.dt);
        title(titlestr)
        legstr = arrayfun(@(x) sprintf('%g', x), ts, 'uniformoutput', false);
        h = legend(legstr, 'location', 'eastoutside');
    end
    
end

function params = setDefaultParameters()
% Set default set of parameters
    
    % Discretization parameters
    params.nx     = 100;
    params.ny     = 10;
    params.dt     = 1e-3;
    params.rampup = 10;
    params.fixedtsteps = [1e-5; 0.01; 0.02; 0.03; 0.04; 0.1; 0.5];
    params.fixedtsteps = [1e-5; 0.01; 0.02; 0.04; 0.08];
    params.totime = params.fixedtsteps(end);
    
    % Flow parameters
    params.perm = 1; % permeability
    params.muW  = 1; % viscosity
    params.poro = 1; % porosity
    params.cW = 0; % compressibility

    % Mechanical parameters
    % Bulk's modulus
    params.K = 1;
    % Poisson's ratio
    params.nu = 0;

end


function plotsol(output)
    
    [ind, xc] = computexc(output);
    ps = output.pressures;
    
    for i = 1 : numel(ps)
        p = ps{i};
        plot(xc, p(ind), 'linewidth', 2);
    end

end

function [ind, xc] = computexc(output)
    
    G = output.G;
    params = output.params;
    ps = output.pressures;

    nx = params.nx;
    ny = params.ny;
    
    ind = (1 : nx)' + floor(ny/2)*nx;
    xc = G.cells.centroids(ind, 1);    
    
end

function plotexactsol(output, dx)
    
    xc = (0 : dx : 1)';
    
    pnorm  = output.pnorm;
    params = output.params;
    cv     = output.cv;
    dtinit = output.dtinit;
    
    ts = params.fixedtsteps;

    nts = numel(ts);
    co = get(gca, 'colororder');

    for i = 1 : nts
        t = ts(i);
        p = pnorm*analyticmandel(cv*t, xc, params, 'num_modes', 1000);
        plot(xc, p, ':', 'color', co(i, :), 'linewidth', 1.4);
    end
end
