function [description, options, state0, model, schedule, plotOptions] = spe10_layer_compositional(varargin)    
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

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
    % One-line description
    description ...
        = ['Water-alternating gas injection in layer 1 of SPE10 Model 2 ', ...
           'using a six-component model. Example from Moncorgé et al, '  , ...
           'J. Comput. Phys., 2018, doi: 10.1016/j.jcp.2018.05.048'      ];
    % Optional input arguments
    options = struct('dt', 25*day); % Target timestep length
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require spe10 ad-core ad-props compositional deckformat
    % Load deck
    fn = fullfile(getDatasetPath('SPE10_Layer_Compositional'), 'TEST1.DATA');
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    % Get grid
    G = SPE10_setup(1);
    G = computeGeometry(G);
    % Get rock
    rock = initEclipseRock(deck);
    rock = compressRock(rock, G.cells.indexMap);
    % Get fluid
    eos  = initDeckEOSModel(deck);
    fluid = initDeckADIFluid(deck);
    inj_comp =  [1.0 0.0 0.0 0.0 0.0 0.0];
    p_std = 14.7*psia;
    T_std = 288.7;
    mc = sum(inj_comp.*eos.CompositionalMixture.molarMass);
    R = 8.3144598;
    rhoL = p_std/(R*T_std);
    rhoL = rhoL.*mc;
    rhoV = rhoL;
    fluid.rhoOS = rhoL;
    fluid.rhoGS = rhoV;
    % Make model
    model = GenericOverallCompositionModel(G, rock, fluid, eos.CompositionalMixture, 'water', true);
    model.AutoDiffBackend = DiagonalAutoDiffBackend();
    % Get schedule
    schedule = convertDeckScheduleToMRST(model, deck);
    nc = numel(schedule.control);
    schedule0 = schedule;
    dt = cell(nc, 1);
    % Set timesteps
    for i = 1:nc
        t_loc = sum(schedule.step.val(schedule.step.control == i));
        dt{i} = rampupTimesteps(t_loc, options.dt, 8);
    end
    schedule.step.val = vertcat(dt{:});
    schedule.step.control = rldecode((1:nc)', cellfun(@numel, dt));
    schedule.s0 = schedule0;
    % Set injection composition
    for i = 1:3
        schedule.control(i).W(1).components = inj_comp;
        schedule.control(i).W(2).components = inj_comp;
    end
    % Get initial state
    z0 = [0.5 0.03 0.07 0.2 0.15 0.05];
    T = 344;
    state0 = initCompositionalState(G, 4000*psia, T, [0, 0, 1], z0, eos);
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', [1,1.83,0.3]  , ...
                   'Projection'        , 'orthographic', ...
                   'View'              , [0, 90]       , ...
                   'Size'              , [500, 800]    };
end