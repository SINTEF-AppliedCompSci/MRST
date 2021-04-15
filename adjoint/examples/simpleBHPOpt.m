%%  Simple Adjoint Test Using BHP-control

mrstModule add adjoint mimetic incomp

% whether or not to show output
verbose = false;
verboseLevel = 0;

%% Define model
nx = 21; ny = 21; nz = 1;
G = cartGrid([nx ny nz], [5*nx 5*ny 1*nz]);
G = computeGeometry(G);

c = G.cells.centroids;
rock.perm  = max(10*sin(c(:,2)/25+.5*cos(c(:,1)/25))-9, .01)*1000*milli*darcy;
rock.poro  = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu' , [1, 5] .* centi*poise, ...
                        'rho', [1014, 859].*kilogram/meter^3, ...
                        'n'  , [2, 2], 'sr', [0, 0], 'kwm', [1, 1]);
fluid  = adjointFluidFields(fluid);


%% Wells and initial rates
radius = .1;
totVol = sum(poreVolume(G, rock));
totTime = 500*day;
W = [];
% Injectors along left side:
nInj = 3; % > 1
pos  = (1 : (ny-1)/(nInj-1) : ny)';
posInj  = round(pos);
for k = 1:nInj
    nm = ['inj', num2str(k)];
    W = addWell(W, G, rock, 1+(posInj(k)-1)*nx, 'Type', 'bhp' , 'Val', 500*barsa, ...
                'Radius', radius, 'Name', nm, 'Comp_i', [1, 0], 'Sign', 1, 'InnerProduct', 'ip_tpf');
end
% Producers along right side:
nProd = 5; % >1
pos  = (1 : (ny-1)/(nProd-1) : ny)';
posProd  = round(pos);
for k = 1:nProd
    nm = ['prod', num2str(k)];
    W = addWell(W, G, rock, nx+(posProd(k)-1)*nx, 'Type', 'bhp' , 'Val', 150*barsa, ...
                'Radius', radius, 'Name', nm, 'Comp_i', [1, 0], 'Sign', -1, 'InnerProduct', 'ip_tpf');
end

%% System components
S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, ...
                     'InnerProduct', 'ip_tpf');
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

%% Initialize
state = initResSol(G, 0.0);
state.wellSol = initWellSol(W, 0);

%% Objective function
objectiveFunction = str2func('simpleNPV');

%% Initialize schedule and controls
numSteps = 10;
schedule = initSchedule(W, 'NumSteps', numSteps, 'TotalTime', ...
                        totTime, 'Verbose', verbose);

% box constraints for each well [min rate, max rate]
box = [repmat([300*barsa 700*barsa], nInj, 1); repmat([100*barsa 200*barsa], nProd, 1)];
controls = initControls(schedule, 'ControllableWells', (1:numel(W)), ...
                                  'MinMax', box, ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', numSteps);

%% Run optimization
[simRes, schedule, controls, out] = optimizeObjective(G, S, W, rock, ...
                                        fluid, state, schedule, ...
                                        controls, objectiveFunction, ...
                                        'gradTol',       1e-3, ...
                                        'objChangeTol',  5e-4, ...
                                        'VerboseLevel', verboseLevel);

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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