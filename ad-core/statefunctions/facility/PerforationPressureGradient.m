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
            gp = gp.dependsOn('bhp', 'wellSol');
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
                % Multiply pressure drop to cut of cross-flowing
                % connections. Note that sign == 0 is ambigious here.
                keep = perfAllowXFlow | (psgn == 1 & vdp <= 0) | (psgn == -1 & vdp >= 0);
                dp = dp.*keep;
            end
        end
    end
end