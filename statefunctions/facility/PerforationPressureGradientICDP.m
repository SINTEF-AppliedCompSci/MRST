classdef PerforationPressureGradientICDP < StateFunction
    
    properties
    end
    
    methods
        
        function gp = PerforationPressureGradientICDP(varargin)
            
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', 'ConnectionPressureDrop'});
            gp = gp.dependsOn('pressure', 'state');
            gp = gp.dependsOn('bhp', 'state');
            gp.label = 'p_c-p_{bh}-g \Delta z \rho_{w}';
            
        end
        
        function dp = evaluateOnDomain(prop, model, state)
            
            [map, cdp] = prop.getEvaluatedDependencies(state, 'FacilityWellMapping', 'ConnectionPressureDrop');
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
            
            cp = bhp + cdp;
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