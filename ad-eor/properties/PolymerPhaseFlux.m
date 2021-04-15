classdef PolymerPhaseFlux < StateFunction
% We use a Todd-Longstaff model. It implies that the mobility of the
% polymer is a non-linear function of the polymer concentration.
            
    properties

    end
    
    methods
        function gp = PolymerPhaseFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('FaceConcentration');
            gp = gp.dependsOn({'FaceMobility', 'PermeabilityPotentialGradient'});
        end
        
        function vP = evaluateOnDomain(prop, model, state)
            [mob, kgrad] = prop.getEvaluatedDependencies(state,...
                'FaceMobility', 'PermeabilityPotentialGradient');
            % compute water phase flux first
            vW     = -mob{1}.*kgrad{1};
            cpf    = prop.getEvaluatedDependencies(state, 'FaceConcentration');
            fluid  = model.fluid;
            mixpar = fluid.mixPar;
            cpbar  = cpf/fluid.cpmax;
            a  = fluid.muWMult(fluid.cpmax).^(1-mixpar);
            vP = vW./(a+(1-a)*cpbar).*cpf;
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
