function varargout = getIncompCapillaryPressure(state, fluid, varargin)
%Undocumented Utility Function

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

    varargout = cell(1, nargout);
    if isfield(fluid, 'properties')
        [varargout{:}] = getPropsLegacy(state, fluid);
    else
        [varargout{:}] = getPropsAD(state, fluid);
    end
end

function [pc, dpc] = getPropsLegacy(state, fluid)
    if isfield(fluid, 'pc')
        if nargout > 1
            [pc, dpc] = fluid.pc(state);
        else
            pc = fluid.pc(state);
        end
    else
        pc = [];
        dpc = [];
    end
end

function [pc, dpc] = getPropsAD(state, fluid)
    if isfield(fluid, 'pcOW')
        getDer = nargout > 1;
        s = state.s(:, 1);
        if getDer
            s = initVariablesAD_diagonal(s);
        end
        pc = fluid.pcOW(s);
        if getDer
            dpc = pc.jac{1}.diagonal; 
        end
    else
        pc = [];
        dpc = [];
    end
end
