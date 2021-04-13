%% Mandel problem
%
% References : 
% @article{verruijt2013theory,
%          title={Theory and problems of poroelasticity},
%          author={Verruijt, Arnold},
%          journal={Delft University of Technology},
%          volume={71},
%          year={2013}
%         } (section 3.2)
%
% and
%
% @article{mikelic2014numerical,
%   title={Numerical convergence study of iterative coupling for coupled flow and geomechanics},
%   author={Mikeli{\'c}, Andro and Wang, Bin and Wheeler, Mary F},
%   journal={Computational Geosciences},
%   volume={18},
%   number={3-4},
%   pages={325--341},
%   year={2014},
%   publisher={Springer}
% }

%% Load required modules
mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

clear all
close all

% Discretization parameters
params.nx     = 100;
params.ny     = 10;
params.dt     = 1e-3;
params.rampup = 1;
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

output = mandelrun(params);

%% plot of the solutions at different pre-selected times

G = output.G;
params = output.params;
ps = output.pressures;

nx = params.nx;
ny = params.ny;

ind = (1 : nx)' + floor(ny/2)*nx;
xc = G.cells.centroids(ind, 1);    

ps = output.pressures;
ts = params.fixedtsteps;

figure 
hold on

for i = 1 : numel(ps)
    p = ps{i};
    plot(xc, p(ind), 'linewidth', 2);
end

% plot the analytical solution (in dotted lines)

dx = 1e-3;
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

legstr = arrayfun(@(x) sprintf('%g', x), ts, 'uniformoutput', false);
h = legend(legstr, 'location', 'eastoutside');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

