classdef FacilityWellMapping < GridProperty
    properties

    end
    
    methods

        function s = evaluateOnDomain(prop, model, state)
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            [wc, p2w] = getActiveWellCells(model, wellSol);
            W = model.getWellStruct(actWellIx);
            isInj = vertcat(W.sign) > 0;
            
            s = struct('active', actWellIx,... % Indices of active wells
                       'cells', wc, ... % Cells where wells are perforated
                       'perf2well', p2w, ... % Perf to well-map
                       'isInjector', isInj, ...
                       'W', W); % Actual well structs
        end
    end
end