function states = addFluxesToRestartStates(states, G, T, rock, fluid, deck)
%Expand Reservoir State to Include Phase Fluxes.

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

% The rock structure is needed when the scaling of relative permeability is
% used.
    
    gravity reset on
    
    % Set up model for property evaluations
    model = selectModelFromDeck(G, rock, fluid, deck);
    % Set up automatically the correct property functions by using validateModel
    model = model.validateModel(); 
    props = model.FlowPropertyFunctions;
    
    % Set up only the operators we need
    model.operators = setupOperatorsTPFA(G, [], 'trans', sparse(numel(T),1), ...
                                            'porv', sparse(G.cells.num,1));
    op  = model.operators;
    gdz = model.getGravityGradient();


    if ~iscell(states)
        states = {states};
    end
    for k =1:numel(states)

        state = states{k};
        state = model.validateState(state);
        [fp, fpname] = props.getPropertyContainer();
        state.(fpname) = fp;
        
        pc  = props.getProperty(model, state, 'CapillaryPressure');  
        rho = props.getProperty(model, state, 'Density');
        mob = props.getProperty(model, state, 'Mobility');
        
        p = model.getProp(state, 'pressure');
        
        p_phase = getPhasePressures(p, pc);
        
        getflux = @(phnum) computeFluxPhase(phnum, rho, p_phase, mob, gdz, T, op);

        % fluxes
        v = cell(1,3);
        phNum = 0;
        if model.water
            phNum = phNum + 1;
            v{phNum} = getflux(phNum);
        end
    
        if model.oil
            phNum = phNum + 1;
            v{phNum} = getflux(phNum);
        end
        
        if model.gas
            phNum = phNum + 1;
            v{phNum} = getflux(phNum);
        end
        
        states{k} = model.storeFluxes(states{k}, v{1}, v{2}, v{3});
    end
end

function v = computeFluxPhase(phNum, rho, p_phase, mob, gdz, T, op)

    rhof  = op.faceAvg(rho{phNum});
    pot   = op.Grad(p_phase{phNum}) - rhof.*gdz;
    upc   = double(pot)<=0;
    dflux = -T.*pot;
    v     = dflux.*op.faceUpstr(upc, mob{phNum});
end
