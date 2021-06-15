classdef WellPhaseFluxTotalFixed < StateFunction
    % Get phase-flux between well-bore and reservoir
    properties
    end
    
    methods
        function gp = WellPhaseFluxTotalFixed(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', 'PressureGradient', 'WellIndex'});
            gp = gp.dependsOn({'Mobility'}, 'FlowPropertyFunctions');
        end
        
        function q_ph = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            mob = model.ReservoirModel.getProps(state, 'Mobility');
            nph = numel(mob);
            q_ph = cell(1, nph);
            W = map.W;
            if isempty(W)
                return
            end
            [dp, wi] = prop.getEvaluatedDependencies(state, 'PressureGradient', 'WellIndex');
            
            ws = state.wellSol(map.active);
            qT = sum(vertcat(ws.flux), 2);
            
            mobw = cellfun(@(x) x(map.cells), mob, 'UniformOutput', false);
            mobt = 0;
            for i = 1:nph
                mobt = mobt + mobw{i};
            end
            isInjector = map.isInjector(map.perf2well);
            
            vTdp = value(-dp);
            injection = vTdp > 0;
            production = ~injection & vTdp ~= 0;
            crossflow = (injection & ~isInjector) | ...
                        (production & isInjector);
            f = cell(1, nph);
            for i = 1:nph
                f{i} = mobw{i}./mobt;
            end
            if any(injection)
                compi = vertcat(W.compi);
                if any(crossflow)
                    dispif(model.verbose > 1, 'Crossflow occuring in %d perforations\n', sum(crossflow));
                    % Compute cross flow for this phase. The approach here
                    % is to calculate (as doubles) the volumetric inflow of
                    % all phases into the well-bore. If a well has
                    % cross-flow, the phase distribution of the
                    % cross-flowing volume is assumed to reflect the inflow
                    % conditions, neglecting density change throughout the
                    % wellbore.
                    q_wb = bsxfun(@times, value(mobw), vTdp);
                    compi = crossFlowMixture(q_wb, compi, map);
                end
                compi_perf = compi(map.perf2well, :);
                for i = 1:nph
                    f{i}(injection) = compi_perf(injection, i);
                end
            end
            for i = 1:nph
                q_ph{i} = f{i}.*qT;
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
