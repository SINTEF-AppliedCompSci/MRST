classdef WellboreComponentPhaseFlux < StateFunction
% State function fore wellbore component phase flux

    methods
        %-----------------------------------------------------------------%
        function cf = WellboreComponentPhaseFlux(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'massFlux'}, 'state');
            cf.label = 'V_{i,\alpha}^w';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function v = evaluateOnDomain(prop, model, state)
            
            % Get mass flux, phase mass per component, and phase upind flag
            vt = state.massFlux;
            [mass, flag] = model.getProps(state, ...
                'ComponentPhaseMass', 'PhaseUpwindFlag');
            nc    = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();
            massT = 0;
            % Compute total mass
            for c = 1:nc
                for ph = 1:nph
                    m = mass{c,ph};
                    if isempty(m), continue; end
                    massT = massT + m;
                end
            end
            
            % Compute phase fluxes as total mass flux multiplied by
            % component mass fractions
            v = cell(nc, nph);
            for c = 1:nc
                for ph = 1:nph
                    m = mass{c,ph};
                    if isempty(m), continue; end
                    % Compute component upstream-weighted fraction of
                    % component c in phase ph
                    x = model.operators.faceUpstr(flag{ph}, m./massT);
                    v{c, ph} = vt.*x;
                end
            end
            
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