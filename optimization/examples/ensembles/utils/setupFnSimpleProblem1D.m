function problem = setupFnSimpleProblem1D(n)
% setup simple problem for testing ensemble optimiation. Adopted from 
% ensemblePackedProblemsExample

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

nx = 100;
rng(nx*n);
G = cartGrid(100, 1000*meter);
G = computeGeometry(G);

% One pore-volumes with maximum porosity
time  = 5*year;
irate = sum(G.cells.volumes)*0.9/time;

poro  = reshape(gaussianField(G.cartDims, [0.25 0.75], [11 3]), [], 1);
perm  = poro.^3.*(1e-5)^2./(0.81*72*(1-poro).^2);
rock = makeRock(G, perm, poro);

dt = rampupTimesteps(time, time/50);
fluid = initSimpleADIFluid('mu', [1, 3]*centi*poise, 'n', [2, 2], 'phases', 'wo', 'cR', 1e-8/barsa, 'c', [1e-6, 1e-5]/barsa);

state0 = initResSol(G, 100*barsa, [0, 1]);

% setup injector and producer
W = addWell([], G, rock,           1, 'Type', 'rate', 'Val', irate,     'Name', 'I', 'comp_i', [1 0], ...
            'Sign',  1, 'lims', struct('bhp', 450*barsa, 'rate', irate));
W = addWell(W,  G, rock, G.cells.num, 'Type', 'bhp',  'Val', 100*barsa, 'Name', 'P', 'comp_i', [0 1], ...
            'Sign', -1, 'lims', struct('lrat', -2*irate, 'bhp', 100*barsa));
schedule = simpleSchedule(dt, 'W', W);

model       = GenericBlackOilModel(G, rock, fluid, 'gas', false);
description = sprintf('Average porosity %1.2f, average perm %1.2f md', mean(poro), mean(perm)/(milli*darcy));
problem     = packSimulationProblem(state0, model, schedule, 'tmp', 'Name', num2str(n), 'Description', description);
end
