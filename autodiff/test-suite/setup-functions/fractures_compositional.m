function setup = fractures_compositional(varargin)
% Setup function for a compositional 2D example with fractures
%
% SYNOPSIS:
%   setup = fractures_compositional('pn1', pv1, ...)
%   setup = fractures_compositional(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Compositional case with fractures, similar to the example from Section
%   5.2 of Moyner & Tchelepi, SPE J, 23(6), doi: 10.2118/182679-PA. The
%   reservoir section consists of five low-permeability bands containing
%   thirteen partially intersecting fractures. Initially, the reservoir
%   contains a mixture consisting of 60% n-decane, 10% carbon dioxide, and
%   30% methane. The reservoir is produced from a well under constant bhp
%   control placed near the northeast corner, supported by a fixed-rate
%   injection of carbon dioxide containing 10% n-decane from a well near
%   the southwest corner.
%
%   Configurable parameters: none
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
    % One-line description
    description ...
        = ['Compositional example with fractures, similar to the ' , ...
           'example from Section 5.2 of Moyner & Tchelepi, SPE J, ', ...
           '23(6), doi: 10.2118/182679-PA'                         ];
       
    % Optional input arguments
    options = struct();
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props compositional
    
    % Name tags for fluid model and case
    fluid_name = 'simple';
    
    % Load the grid and set up petrophysical model
    pth  = fullfile(getDatasetPath('MSFractures'), 'setup_fracture.mat');
    data = load(pth);
    G    = data.G;
    perm = data.perm;
    rdim = [1000 500];
    G.nodes.coords = bsxfun(@times, G.nodes.coords, rdim);
    G    = computeGeometry(G);
    rock = makeRock(G, perm*milli*darcy, 0.3);
    rock.perm(G.cells.tag > 0) = 10*darcy;
    pv   = poreVolume(G, rock);
    
    % Define fluid and PVT model
    [cf, info] = getBenchmarkMixture(fluid_name);
    eos   =  EquationOfStateModel(G, cf);
    fluid = initSimpleADIFluid('rho', [1000, 500, 500]     , ...
                               'mu' , [1, 1, 1]*centi*poise, ...
                               'n'  , [2, 2, 2]            , ...
                               'c'  , [1e-5, 0, 0]/barsa   );
    
    % Setup the wells
    minP    = 50*barsa;
    resP    = info.pressure;
    totTime = 7*year;
    irate   = 300*0.25*sum(pv)/totTime;
    icell   = findEnclosingCell(G, [0.025, 0.05].*rdim);
    pcell   = findEnclosingCell(G, [0.975, 0.95].*rdim);
    W = addWell([], G, rock, icell, 'comp_i', [0, 1], ...
                'name', 'I', 'Type', 'rate', 'Val', irate);
    W = addWell(W, G, rock, pcell, 'comp_i', [0.5, 0.5], ...
                'Name', 'P', 'Val', minP, 'sign', -1);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    
    % Build the model and set accelerated AD backend
    model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
    
    % Set the initial state, which must be expanded with information about
    % components and temperature
    state0   = initResSol(G, resP, [1 0]);
    state0.T = repmat(info.temp, G.cells.num, 1);    
    state0.components = repmat(info.initial, G.cells.num, 1);
    
    % Build the simulation schedule: we run with uniform time steps of 20
    % days, which give a reasonable compromise between having too high CFL
    % numbers in the fractures and too low CFL numbers in the background
    % matrix. In addition, we add a standard rampup to stabilize the
    % displacement fronts as they move into the reservoir.
    dt       = rampupTimesteps(totTime, 20*day);
    schedule = simpleSchedule(dt, 'W', W);
    
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', [2,1,1], 'Size', [800, 500]};
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end
