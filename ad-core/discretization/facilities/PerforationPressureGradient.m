classdef PerforationPressureGradient < AutoDiffFunction
    properties

    end
    
    methods
        function gp = PerforationPressureGradient(varargin)
            gp@AutoDiffFunction(varargin{:});
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
        end
    end
end