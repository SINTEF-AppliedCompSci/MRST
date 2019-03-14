classdef PerforationPressureGradient < GridProperty
    properties

    end
    
    methods
        function gp = PerforationPressureGradient(varargin)
            gp@GridProperty(varargin{:});
        end
        function dp = evaluateOnDomain(prop, model, state)
            map = model.getProp(state, 'FacilityWellMapping');
            wellSol = state.wellSol;
            p = model.ReservoirModel.getProps(state, 'pressure');
            
            % Temporary
            isbhp = strcmpi(state.FacilityState.names, 'bhp');
            bhp = state.FacilityState.primaryVariables{isbhp}(map.perf2well);
            cp = bhp + vertcat(wellSol(map.active).cdp);
            dp = p(map.cells) - cp;
        end
    end
end