classdef PerforationPressureGradient < StateFunction
    % Calculate the pressure difference between the reservoir cells and the
    % well-bore pressure (assuming pre-computed static connection
    % pressure-drop in wellSol.cdp)
    properties

    end
    
    methods
        function gp = PerforationPressureGradient(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FacilityWellMapping');
            gp = gp.dependsOn('pressure', 'state');
            gp = gp.dependsOn('bhp', 'state');
            gp.label = 'p_c-p_{bh}-g \Delta z \rho_{w}';
        end
        function dp = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            wellSol = state.wellSol;
            p = model.ReservoirModel.getProps(state, 'pressure');
            
            % Temporary
            isbhp = strcmpi(state.FacilityState.names, 'bhp');
            if any(isbhp)
                bhp = state.FacilityState.primaryVariables{isbhp}(map.perf2well);
            else
                bhp = vertcat(wellSol(map.active).bhp);
                bhp = bhp(map.perf2well);
            end
            cp = bhp + vertcat(wellSol(map.active).cdp);
            dp = p(map.cells) - cp;
            allowXFlow = cellfun(@(x) x.allowCrossflow, model.WellModels(map.active));
            if ~all(allowXFlow)
                % We have cross-flow, potentially
                vdp = value(dp);
                sgn = vertcat(map.W.sign);
                perfAllowXFlow = allowXFlow(map.perf2well);
                psgn = sgn(map.perf2well);
                % Set pressure drop to zero to cut of cross-flowing
                % connections. Note that sign == 0 is ambigious here.
                keep = perfAllowXFlow | (psgn == 1 & vdp <= 0) | (psgn == -1 & vdp >= 0);
                if isa(dp, 'ADI')
                    dp.val(~keep) = 0;
                else
                    dp(~keep) = 0;
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
