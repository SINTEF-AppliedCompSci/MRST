function [state0, model, schedule, ref] = setupSimpleCompositionalExample(useNatural)
%Undocumented Utility Function

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

    if nargin == 0
        useNatural = true;
    end
    require deckformat compositional ad-props
    pth = getDatasetPath('simplecomp');
    fn  = fullfile(pth, 'SIMPLE_COMP.DATA');
    % Read deck
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    % Set up grid
    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    % Set up rock
    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    fluid = initDeckADIFluid(deck);
    % Define some surface densities
    fluid.rhoOS = 800;
    fluid.rhoGS = 10;

    eos = initDeckEOSModel(deck);
    arg = {G, rock, fluid, eos, 'water', false};
    if useNatural
        model = GenericNaturalVariablesModel(arg{:});
    else
        model = GenericOverallCompositionModel(arg{:});
    end
    schedule = convertDeckScheduleToMRST(model, deck);

    % Manually set the injection composition
    [schedule.control.W.components] = deal([0, 1, 0]);
    % Injection is pure gas
    [schedule.control.W.compi] = deal([1, 0]);

    % Set up initial state
    % The problem is defined at 150 degrees celsius with 75 bar initial
    % pressure. We set up the initial problem and make a call to the flash
    % routines to get correct initial composition.

    for i = 1:numel(schedule.control.W)
        schedule.control.W(i).lims = [];
    end

    % Initial conditions
    z0 = [0.6, 0.1, 0.3];
    T = 150 + 273.15;
    p = 75*barsa;

    state0 = initCompositionalState(G, p, T, [1, 0], z0, eos);
    if nargout > 3
        ref = load(fullfile(pth, 'comparison.mat'));
    end
end
