classdef PolymerPhaseFlux < StateFunction
    properties

    end
    
    methods
        function gp = PolymerPhaseFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseFlux', 'FaceConcentration'});
        end
        
        function v = evaluateOnDomain(prop, model, state)
            [phaseFlux, cpf] = prop.getEvaluatedDependencies(state,...
                'PhaseFlux', 'FaceConcentration');
            nph = numel(phaseFlux);
            v = cell(1, nph+1);
            for i = 1:nph
                v{i} = deal(phaseFlux{i});
            end
            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cpbar   = cpf/fluid.cpmax;
            a = fluid.muWMult(fluid.cpmax).^(1-mixpar);
            v{nph+1} = v{1}./(1+(1-a)*cpbar).*cpf;
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
