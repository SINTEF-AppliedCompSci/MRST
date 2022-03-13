function [sol, report] = migrateInjection(Gt, traps, petrodata, wellCell, varargin)
% Run a simple injection scenario and visualize each time step
%
% SYNOPSIS:
%   function sol = migrateInjection(Gt, traps, petrodata, wellCell, varargin)
%
% DESCRIPTION:
%
% This script runs a simple injection scenario based on a top-surface grid
% and a single injector well located in a specified grid cell.  The
% simulation uses an incompressible fluid with linear relative permeabilities
% and sharp interface.  Rock is incompressible and homogeneous, with
% permeability and porosity provided by the 'petrodata' argument.  After
% computation of a timestep, the variously trapped volumes are computed, and
% the result is visualized (either using 'plotPanelVE', or simply by using
% 'plotCellData').  
%
% For the simulation, pressure is solved using two-point flux approximation.
% Saturations are computed using implicit transport.
%
% PARAMETERS:
%   Gt        - Top surface grid 
%   traps     - trapping structure object, as returned by a call to
%               'findTrappingStructure'.  Only used for visualizations
%               involving plotPanelVE.  Can be left empty, in which case it
%               will be computed internally if needed.
%   petrodata - structure with the fields 'avgperm' (average permeability)
%               and 'avgporo' (average porosity).  These are single scalars,
%               which will be attributed uniformly to all cells.
%   wellCell  - Index of the cell containing the injection well.
%   varargin  - Allows special options to be set / defaults to be
%               overridden.  For instance, the default injection and
%               migration times can be changed, as well as the number of timesteps.
%
% RETURNS:
%   sol - simulation solution for the last timestep 
%   report - struct reporting the CPU time for the splitting steps
%
% SEE ALSO:
%

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
require incomp
%%% Process options
opt = struct('amount',      1, ...
            'T_injection',  100*year,  ...
            'T_migration',  1000*year, ...
            'topPressure',  300*barsa, ...
            'Ni',           20,        ...
            'Nm',           40,        ...
            'plotPanel',    false,     ...
            'view',         [-85 70],...
            'plotPlume',    true, ...
            'plotHist',     false,     ...
            'wireH',        true,...
            'wireS',        true);
opt = merge_options(opt, varargin{:});

%%% Precompute traps if missing and required
if opt.plotPanel && isempty(traps)
  traps = trapAnalysis(Gt, false);
end

%%% Take care of global variables
% When launched from the interactive viewer, we use a global variable to
% indicate when the user has pushed the 'Abort' button. If launched as a
% standalone simulation, this variable must be set to false.
global veSimAborted;
if isempty(veSimAborted)
   veSimAborted = false;
end

%%% Set up fluid
muw = 0.30860;  rhow = 975.86; sw    = 0.1;
muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
kwm = [0.2142 0.85];

mu  = [muc  muw ] .* centi*poise;
rho = [rhoc rhow] .* kilogram/meter^3;

fluid = initSimpleVEFluid_s('mu' , mu , 'rho', rho, ...
                            'height'  , Gt.cells.H, ...
                            'sr', [srco2, sw],'kwm',kwm);
fluid.res_gas = srco2;
fluid.res_water = sw;

%%% Define vertical pressure distribution
gravity on;
grav     = gravity();
topPos   = min(Gt.cells.z);
pressure = @(z) opt.topPressure + rho(2)*(z - topPos)*grav(3);

%%% Schedule
T_tot = opt.T_injection + opt.T_migration;
dTi   = opt.T_injection / opt.Ni;  % short time steps during injection
dTm   = opt.T_migration / opt.Nm;  % longer steps during migration


%%% Set up rock properties and compute transmissibilities
% We use the averaged values for porosity and permeability as given in the
% Atlas tables. Since cellwise data is not present, we assume to averaged
% values to be valid everywhere.
G = Gt.parent;
pd = petrodata;
rock.poro = repmat(pd.avgporo, G.cells.num, 1);
rock.perm = repmat(pd.avgperm, G.cells.num, 1);
rock2D    = averageRock(rock, Gt);
T = computeTrans(Gt, rock2D);
T = T.*Gt.cells.H(gridCellNo(Gt));

%%% Set up well and boundary conditions
% This example is using an incompressible model for both rock and fluid. If
% we assume no-flow on the boundary, this will result in zero flow from a
% single injection well. However, this can be compensated if we use the
% physical understanding of the problem to set appropriate boundary
% conditions: The atlas formations are enormous compared to the volume of
% the injected CO2. Thus, it is impossible that the injection will change
% the composition of the formation significantly. We therefore assume that
% the boundary conditions can be set equal to hydrostatic pressure to drive
% flow.

% Add an injector well for the CO2 with Mt annual injection
rate = opt.amount*1e9*kilogram/(year*rhoc*kilogram*meter^3);
W = addWell([], G, rock, wellCell,...
   'Type', 'rate', 'Val', rate, 'comp_i', [1,0], 'name', 'Injector', ...
   'InnerProduct', 'ip_tpf');

% Add pressure boundary
bnd = boundaryFaces(Gt);
bc = addBC([], bnd, 'pressure', pressure(Gt.faces.z(bnd)), 'sat', [0 1]);

% Convert to 2D wells
W2D = convertwellsVE(W, G, Gt, rock2D,'ip_tpf');

%%%  Set up initial reservoir conditions
% The initial pressure is set to hydrostatic pressure. Setup and plot.
sol = initResSolVE_s(Gt, pressure(Gt.cells.z), 0);
sol.wellSol = initWellSol(W2D, 0);
sol.h = zeros(Gt.cells.num, 1);

[ii, jj] = ind2sub(G.cartDims, G.cells.indexMap);

opts = {'slice',     double([ii(wellCell), jj(wellCell)]),...
        'Saxis',     [0 1-fluid.res_water], ...
        'view',      opt.view,...
        'plotPlume', opt.plotPlume, ...
        'wireH',     opt.wireH,...
        'wireS',     opt.wireS,...
        'maxH',      (max(Gt.cells.z) - min(Gt.cells.z))/3, ...
        'plotHist',  opt.plotHist};

if opt.plotPanel
    fh = plotPanelVE(G, Gt, W, sol, 0, ...
       [volumesVE(Gt, sol, rock2D, fluid, traps) 0], opts{:});
else
    fh = figure;
end

%%% Run the simulation
% Solve for all timesteps, and plot the results at each timestep.
t = 0;
tt = ' (Injecting)';
totVol = 0;

% Estimate total timesteps
tstep_tot = ceil(opt.T_injection/dTi) + ceil(opt.T_migration/dTm);
i = 1;

% Waitbar to show progress
hwbar = waitbar(0, 'Initializing...');
wbar = @(i, t, status) waitbar(t/T_tot, hwbar, ...
   sprintf('Timestep %d of %d, T=%s%s', i, tstep_tot, ...
   formatTimeRange(floor(t/year)*year), status));
[ctime, cputimeT,cputimeP] = deal(0);
while t < T_tot
    if ishandle(hwbar)
        wbar(i, t, tt);
    end
    if t >= opt.T_injection
        W2D = [];
        dT  = dTm;
        % Do a shorter timestep
        dT = min(dT, T_tot - t);
        tt  = ' (Migrating)';
    else
        dT = dTi;
        dT = min(dT, opt.T_injection - t);
        W2D(1).cells = double(wellCell);
    end
    
    tStartP = tic;
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc);
    cputimeP(i) = toc(tStartP);
    tStartT = tic;
    sol = implicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false, 'Trans', T);
    cputimeT(i) = toc(tStartT);
    t = t + dT;
    ctime(i)=t;

    % Plotting
    [s, h, hm] = normalizeValuesVE(Gt, sol, fluid);
    sol.h = h;
    sol.h_max = hm;
    if ~ishandle(fh) || veSimAborted
        disp 'Simulation aborted!'
        break
    end

    if opt.plotPanel
        % Use advanced plotting
        if ~isempty(W2D)
            totVol = totVol + W2D.val*dT;
        end
        vol = volumesVE(Gt, sol, rock2D, fluid, traps);
        plotPanelVE(G, Gt, W, sol, t, [vol totVol], opts{:});

    else
        set(0,'CurrentFigure', fh);
        [a,b] = view();
        clf
        plotCellData(G, s, 'edgec', 'k', 'edgea', .1, 'edgec', [.6 .6 .6]);
        plotWell(G, W); caxis([0 .9]);
        title([formatTimeRange(t) tt])
        colorbar
        axis tight off
        view(a,b);
        drawnow
    end
    i = i + 1;
end
fprintf('Total CPU time: %f %f (%d steps)\n',sum(cputimeP),sum(cputimeT), i-1);
report.transport = cputimeT;
report.pressure = cputimeP;
report.time = ctime;
if ishandle(hwbar), close(hwbar); end
end
