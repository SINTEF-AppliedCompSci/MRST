function  [model, initState, schedule] = setupAquifertest(fn)
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

    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    % G = initEclipseGrid(deck, 'SplitDisconnected', false);
    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    % rock.perm = 1e8*rock.perm;

    % end-point scaling, so include grid
    fluid = initDeckADIFluid(deck, 'G', G);

    output = processAquifer(deck, G);
    aquifers     = output.aquifers;
    aquind       = output.aquind; 
    initval      = output.initval;
    aquiferprops = output.aquiferprops;

    model = AquiferBlackOilModel(G, rock, fluid, aquifers, aquind, aquiferprops, ...
                                 'modeltype', 'oilwater');

    if isfield(deck.RUNSPEC, 'VAPOIL')
        model.vapoil = deck.RUNSPEC.VAPOIL;
        model.vapoil = deck.RUNSPEC.VAPOIL;
    else
        model.vapoil = 0;        
        model.vapoil = 0;        
    end
    if isfield(deck.RUNSPEC, 'DISGAS')
        model.disgas = deck.RUNSPEC.DISGAS;
        model.disgas = deck.RUNSPEC.DISGAS;
    else
        model.disgas = 0;        
        model.disgas = 0;        
    end
    model.FacilityModel = selectFacilityFromDeck(deck, model);
    % nc = G.cells.num;
    % initState.pressure = 250*barsa*ones(nc, 1);
    % s = 0.3*ones(nc, 2);
    % s(:, 2) = 1 - s(:, 1);
    % initState.s = s;
    initState = initStateDeck(model, deck);
    initState.aquiferpressures = initval.pressures;
    initState.aquifervolumes   = initval.volumes;

    schedule = convertDeckScheduleToMRST(model, deck);

end
