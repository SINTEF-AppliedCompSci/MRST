classdef FacilityWellMapping < GridProperty
    properties

    end
    
    methods

        function s = evaluateOnDomain(prop, model, state)
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            wc = getActiveWellCells(model, wellSol);
            W = model.getWellStruct(actWellIx);
            p2w = getPerforationToWellMapping(W);
            isInj = vertcat(W.sign) > 0;
            
            wsum = sparse(p2w, (1:numel(p2w))', 1);
            
            s = struct('active', actWellIx,... % Indices of active wells
                       'cells', wc, ... % Cells where wells are perforated
                       'perf2well', p2w, ... % Perf to well-map
                       'isInjector', isInj, ...
                       'perforationSum', wsum, ...
                       'W', W); % Actual well structs
        end
    end
end