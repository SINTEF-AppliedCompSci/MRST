classdef PhaseThermalEnergy < StateFunction
%State function for thermal energy in each fluid phase
    
    properties
        useEnthalpy = false;
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseThermalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ComponentPhaseMass'});
            if gp.useEnthalpy
                gp = gp.dependsOn({'PhaseEnthalpy'}, 'PVTPropertyFunctions');
            else
                gp = gp.dependsOn({'PhaseInternalEnergy'});
            end
            gp.label = 'E_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function energy = evaluateOnDomain(prop,model, state)
            mass = prop.getEvaluatedDependencies(state, 'ComponentPhaseMass');
            if prop.useEnthalpy
                e = model.getProp(state, 'PhaseEnthalpy');
            else
                e = prop.getEvaluatedDependencies(state, 'PhaseInternalEnergy');
            end
            nph    = model.getNumberOfPhases();
            ncomp  = model.getNumberOfComponents();
            energy = cell(1, nph);
            for i = 1:nph
                phmass = 0;
                for j = 1:ncomp
                    if isempty(mass{j,i}), continue; end
                    phmass = phmass + mass{j,i};
                end
                energy{i} = phmass.*e{i};
            end
        end

    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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