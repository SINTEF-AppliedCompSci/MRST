function [deck, schedule, model, rock] = setupADcase(fn)
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

    moddir = mrstPath('query', 'ad-testdata');
    fn = fullfile(moddir, fn);
    if ~exist(fn, 'file')
        error(['Did not find dataset at expected location: (', fn , ')'])
    end

    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);

    fluid = initDeckADIFluid(deck);
    
    model = selectModelFromDeck(G, rock, fluid, deck);
    schedule = convertDeckScheduleToMRST(model, deck);
    
end
