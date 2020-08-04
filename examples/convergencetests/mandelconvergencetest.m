function output = mandelconvergencetest
    
    params = setDefaultParameters();
    
    ts = params.fixedtsteps;
    
    % Spatial convergence
    nxs = [5];
    nnxs = numel(nxs);
    states = cell(nnxs, 1);
    
    close all
    
    for i = 1 : nnxs
        params.nx = nxs(i);
        output{i} = mandelrun(params);
        figure
        hold on
        plotsol(output{i}, params);
        plotexactsol(output{i}, params, 1e-2);
        titlestr = sprintf('spatial convergence, nx = %d', params.nx);
        title(titlestr)
        legstr = arrayfun(@(x) sprintf('%g', x), ts, 'uniformoutput', false);
        h = legend(legstr, 'location', 'eastoutside');
    end
    
    output = [];
    
end

function params = setDefaultParameters()
% Set default set of parameters
    
    % Discretization parameters
    params.nx     = 100;
    params.ny     = 10;
    params.dt     = 1e-3;
    params.rampup = 10;
    params.fixedtsteps = [1e-5; 0.01; 0.02; 0.03; 0.04; 0.1; 0.5];
    params.fixedtsteps = [1e-5; 0.01; 0.02; 0.03; 0.04];
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


function plotsol(output, params)
    
    [ind, xc] = computexc(output, params);
    ps = output.pressures;
    
    for i = 1 : numel(ps)
        p = ps{i};
        plot(xc, p(ind), 'linewidth', 2);
    end

end

function [ind, xc] = computexc(output, params)
    
    nx = params.nx;
    ny = params.ny;
    G = output.G;
    ps = output.pressures;
    
    ind = (1 : nx)' + floor(ny/2)*nx;
    xc = G.cells.centroids(ind, 1);    
    
end

function plotexactsol(output, params, dx)
    
    xc = (0 : dx : 1)';
    ts = params.fixedtsteps;
    
    pnorm  = output.pnorm;
    cv     = output.cv;
    dtinit = output.dtinit;
    
    nts = numel(ts);
    co = get(gca, 'colororder');

    for i = 1 : nts
        t = ts(i);
        p = pnorm*analyticmandel(cv*t, xc, params, 'num_modes', 1000);
        plot(xc, p, ':', 'color', co(i, :), 'linewidth', 1.4);
    end
end
