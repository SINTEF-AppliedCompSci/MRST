classdef FacilityWellMapping < StateFunction
    properties

    end
    
    methods
        function gp = FacilityWellMapping(varargin)
            gp@StateFunction(varargin{:});
        end
        function s = evaluateOnDomain(prop, model, state)
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            wc = getActiveWellCells(model, wellSol);
            W = model.getWellStruct(actWellIx);
            p2w = getPerforationToWellMapping(W);
            if isempty(W)
                isInj = [];
            else
                isInj = vertcat(W.sign) > 0;
            end
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