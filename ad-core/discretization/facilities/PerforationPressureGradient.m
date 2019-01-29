classdef PerforationPressureGradient < GridProperty
    properties

    end
    
    methods

        function dp = evaluateOnDomain(prop, model, state)
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            [wc, p2w] = getActiveWellCells(model, wellSol);
            
            p = model.ReservoirModel.getProps(state, 'pressure');
            
            % Temporary
            isbhp = strcmpi(state.FacilityState.names, 'bhp');
            bhp = state.FacilityState.primaryVariables{isbhp}(p2w);
            cp = bhp + vertcat(wellSol(actWellIx).cdp);
            dp = p(wc) - cp;
        end
    end
end