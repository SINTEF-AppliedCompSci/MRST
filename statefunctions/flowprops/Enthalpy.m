classdef Enthalpy < StateFunction
%State function for enthaply

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = Enthalpy(model, varargin)
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'enthalpy'
                    % Enthalpy is the primary variable
                case 'temperature'
                    % Temperature is the primary variable
                    gp = gp.dependsOn({'PhasePressures', 'Density'});
                    gp = gp.dependsOn('PhaseInternalEnergy', 'FlowPropertyFunctions');
            end
            gp.label = 'h_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function h = evaluateOnDomain(prop, model, state)
            switch model.thermalFormulation
                case 'enthalpy'
                    error('This line should not be reached - something is wrong');
                case 'temperature'
                    [p, rho] = prop.getEvaluatedDependencies(state, 'PhasePressures'     , ...
                                                                    'Density'            );
                    u = model.getProps(state, 'PhaseInternalEnergy');
                    nph = model.getNumberOfPhases();
                    h = cell(1, nph);
                    for i = 1:nph
                        h{i} = u{i} + p{i}./rho{i};
                    end
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