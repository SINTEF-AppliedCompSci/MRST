classdef TotalThermalEnergy < ComponentTotalMass
%State function for computing the total thermal energy in each cell
    
    methods
        function gp = TotalThermalEnergy(model, varargin)
            gp = gp@ComponentTotalMass(model, varargin{:});
            gp.dependencies = {};
            gp = gp.dependsOn({'PhaseThermalEnergy', 'RockThermalEnergy'});
            gp.label = 'E';
        end
        
        function energy = evaluateOnDomain(prop, model, state)
            [phaseEnergy, rockEnergy] ...
                = prop.getEvaluatedDependencies(state, 'PhaseThermalEnergy', ...
                                                       'RockThermalEnergy' );
            energy = rockEnergy;
            nph = model.getNumberOfPhases();
            for i = 1:nph
                energy = energy + phaseEnergy{i};
            end
            energy = prop.ensureMinimumDerivatives({energy});
            energy = energy{1};
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