classdef UpwindFunctionWrapperDiscretization < UpwindDiscretization
    % Simple wrapper for function handles for upwinding (classical MRST
    % style)
    properties (Access = protected)
        function_handle
    end
    
    methods
        function ufn = UpwindFunctionWrapperDiscretization(model)
            ufn@UpwindDiscretization(model);
            if isempty(model.operators)
                N = getNeighbourship(model.G);
                nc = model.G.cells.num;
                nf = size(N, 1);
                up = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);
            else
                up = model.operators.faceUpstr;
            end
            assert(isa(up, 'function_handle'));
            ufn.function_handle = up;
        end
        
        function v = faceUpstream(wrapper, model, state, flag, cellvalue)
            v = wrapper.function_handle(flag, cellvalue);
        end
        
        function ufn = setFunctionHandle(ufn, up)
            ufn.function_handle = up;
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
