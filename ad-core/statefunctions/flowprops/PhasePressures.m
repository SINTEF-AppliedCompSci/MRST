classdef PhasePressures < StateFunction
    % Pressure for each phase. Will always return one value for each phase,
    % even if all phases have the same pressure.
    properties
        has_pc
    end
    
    methods
        function gp = PhasePressures(model, varargin)
            gp@StateFunction(model, varargin{:});
            f = model.fluid;
            gp.has_pc = (model.water && model.oil && isfield(f, 'pcOW')) || ...
                        (model.gas   && model.oil && isfield(f, 'pcOG')) || ...
                        (~model.oil  &&              isfield(f, 'pcWG'));
            if gp.has_pc
                gp = gp.dependsOn({'CapillaryPressure'}, 'FlowPropertyFunctions'); 
            end
            gp = gp.dependsOn({'pressure'}, 'state');                   
            gp.label = 'p_\alpha';
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            p = model.getProps(state, 'Pressure');
            nph = model.getNumberOfPhases();
            p_phase = cell(1, nph);
            [p_phase{:}] = deal(p);
            if prop.has_pc
                pc = prop.getEvaluatedExternals(model, state, 'CapillaryPressure');
                for i = 1:nph
                    if ~isempty(pc{i})
                        p_phase{i} = p_phase{i} + pc{i};
                    end
                end
            end
        end
    end
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
