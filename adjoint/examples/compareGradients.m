%% Compare Gradient Computed Numerically and by Adjoint Equations
mrstModule add adjoint mimetic incomp

%% Setup model
verbose = true;
verboseLevel = 1;

nx = 10; ny = 10; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100*milli*darcy;
rock.perm = ones( size(rock.perm) )*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu' , [   1,   5].*centi*poise    , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]                 , ...
                        'sr' , [   0,   0]                 , ...
                        'kwm', [   1,   1]);

fluid  = adjointFluidFields(fluid);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose);

%% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

%% Introduce wells
radius = .1;
W = addWell([], G, rock, 1         , 'Type', 'bhp' , 'Val',  100*barsa, 'Radius', radius, 'Name', 'i1', 'Comp_i', [1, 0], 'InnerProduct', 'ip_tpf');
W = addWell( W, G, rock, nx        , 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p1', 'Comp_i', [0, 1], 'InnerProduct', 'ip_tpf');
W = addWell( W, G, rock, nx*ny-nx+1, 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p3', 'Comp_i', [0, 1], 'InnerProduct', 'ip_tpf');
W = addWell( W, G, rock, nx*ny     , 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p3', 'Comp_i', [0, 1], 'InnerProduct', 'ip_tpf');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

resSolInit = initResSol(G, 0.0);
resSolInit.wellSol = initWellSol(W, 0);


totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 10, 'TotalTime', totVol*day, 'Verbose', verbose);


controls = initControls(schedule, 'ControllableWells', (2:4), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 10);


%% Forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

%% Adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
grad   = computeGradient(W, adjRes, schedule, controls);

numGrad = computeNumericalGradient(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction)
adjGrad = cell2mat(grad)


%% Plot results
figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')


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

