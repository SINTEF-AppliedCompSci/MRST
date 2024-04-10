classdef WellboreSurfaceRate < StateFunction
%State function for wellbore surface rate.

    methods
        %-----------------------------------------------------------------%
        function cf = WellboreSurfaceRate(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'ComponentPhaseFlux', 'PhaseUpwdinFlag'}, 'FlowPropertyFunctions');
            cf = cf.dependsOn({'SurfaceDensity'}, 'PVTpropertyFunctions');
            cf.label = 'q_s';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function qs = evaluateOnDomain(prop, model, state)
            
            [q, rhoS, flag] = model.getProps(state, ...
                'ComponentPhaseFlux', 'SurfaceDensity', 'PhaseUpwindFlag');
            [~, iif] = model.getInletSegments();
            rhoS = cellfun(@(flag, rho) ...
                model.parentModel.operators.faceUpstr(flag, rho), ...
                flag, rhoS, 'UniformOutput', false);
            qs = cellfun(@(q, rho) q(iif)./rho(iif), ...
                    q, rhoS, 'UniformOutput', false);
            qs = cellfun(@(qs) model.parentModel.operators.fluxSum(qs), ...
                    qs, 'UniformOutput', false);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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