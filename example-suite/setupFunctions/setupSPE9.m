function [G, rock, fluid, deck, state] = setupSPE9()
%Setup the SPE9 benchmark case for simulation
%
% SYNOPSIS
%   [G, rock, fluid, deck, state] = setupSPE9()
%
% PARAMETERS:
%   none
%
% RETURNS:
%   G     - Grid data structure
%   rock  - Petrophysical variables
%   fluid - Fluid structure
%   deck  - Keywords from the input data file
%   state - Initial reservoir state
%
% SEE ALSO:
%   setupSPE1

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
    mrstModule add deckformat ad-props
    % Read and process file.
    pth = getDatasetPath('spe9');
    fn  = fullfile(pth, 'BENCH_SPE9.DATA');

    deck = readEclipseDeck(fn);

    % The deck is given in field units, MRST uses metric.
    deck = convertDeckUnits(deck);

    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);

    % Create a special ADI fluid which can produce differentiated fluid
    % properties.
    fluid = initDeckADIFluid(deck);

    % The case includes gravity
    gravity reset on

    p0  = deck.SOLUTION.PRESSURE;
    sw0 = deck.SOLUTION.SWAT;
    sg0 = deck.SOLUTION.SGAS;
    s0  = [sw0, 1-sw0-sg0, sg0];
    % Gas in oil phase
    if isfield(deck.SOLUTION, 'RS')
        rs0 = deck.SOLUTION.RS;
    else
        rs0 = 0;
    end
    % Oil in gas phase
    if isfield(deck.SOLUTION, 'RV')
        rv0 = deck.SOLUTION.RV;
    else
        rv0 = 0;
    end
    
    state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
end