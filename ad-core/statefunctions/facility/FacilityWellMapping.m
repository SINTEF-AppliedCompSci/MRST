classdef FacilityWellMapping < StateFunction
    % Generate a struct containing a bunch of useful mappings for a given
    % set of active wells
    properties

    end
    
    methods
        function gp = FacilityWellMapping(varargin)
            gp@StateFunction(varargin{:});
            gp.label = 'W\rightarrow c';
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
