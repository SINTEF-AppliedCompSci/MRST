function [schedule, model, state0] = getBenchmarkAD(name)
    mrstModule add deckformat ad-fi ad-core ad-blackoil ad-props
    name = lower(name);
    gravity reset on
    switch name
        case 'spe1'
            % SPE1 - black oil with gas dissolved in oil.
            [schedule, model, state0] = readAndSetup('SPE', 'SPE1', 'BENCH_SPE1.DATA');
        case 'spe3'
            % SPE3 - black oil with oil dissolved in gas phase.
            [schedule, model, state0] = readAndSetup('SPE', 'SPE3', 'BENCH_SPE3.DATA');
        case 'spe9'
            % SPE9 - black oil with gas dissolaved in oil phase and
            % multiple wells with changing controls.
            [schedule, model, state0] = readAndSetup('SPE', 'SPE9', 'BENCH_SPE9.DATA');
            model.drsMaxRel = .2;
            model.dpMaxRel  = .2;
            model.dsMaxAbs  = .05;
            
        case 'egg'
            % Egg benchmark model from TU Delft (single realization). Uses
            % nearly incompressible oil/water system.
            [schedule, model] = readAndSetup('external', 'TUDelft-EGG', 'BENCH_EGG.DATA');
            G = model.G;
            fluid = model.fluid;
            % Approximate initial conds:
            pr   = 400*barsa;
            rz   = G.cells.centroids(1,3);
            dz   = G.cells.centroids(:,3) - rz;
            rhoO    = fluid.bO(400*barsa)*fluid.rhoOS;
            rhoW    = fluid.bW(400*barsa)*fluid.rhoWS;
            rhoMix  = .1*rhoW + .9*rhoO;
            p0   = pr + norm(gravity)*rhoMix*dz;      
            
            state0 = initResSol(G, p0, [0.1, .90]);            
        case 'simpleow'
            % Simple oil/water test.
            [schedule, model] = readAndSetup('SINTEF', 'simpleOW', 'simple10x1x10.data');
            state0 = initResSol(model.G, model.inputdata.PROPS.PVCDO(1), [.15, .85]);
            
        case 'simplepolymer'
            % Simple polymer test
            [schedule, model] = readAndSetup('SINTEF', 'simplePolymer', 'POLYMER.DATA');
            G = model.G;
            ijk = gridLogicalIndices(G);
            state0 = initResSol(G, model.inputdata.PROPS.PVCDO(1), [ .9, .1]);
            state0.s(ijk{3} == 1, 2) = .9;
            state0.s(ijk{3} == 2, 2) = .8;
            state0.s(:,1) = 1 - state0.s(:,2);
            
            % Add zero polymer concentration to the state.
            state0.c    = zeros(G.cells.num, 1);
            state0.cmax = zeros(G.cells.num, 1);
        otherwise
            error(['Unknown benchmark case ''', name, '''.'])
    end
end

function [schedule, model, state0] = readAndSetup(varargin)
    fn = fullfile(varargin{:});
    mrstModule add ad-unittest
    [deck, schedule, model] = setupADcase(fn);
    if nargout > 2
        state0 = stateFromDeckSolution(deck);
    end
end

function state0 = stateFromDeckSolution(deck)
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
    
    state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
end

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
