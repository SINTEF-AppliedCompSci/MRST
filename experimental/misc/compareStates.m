function stateDiff = compareStates(state1, state2, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('relative', false);
    opt = merge_options(opt, varargin{:});

    stateDiff = state1;
    flds = fieldnames(state1);
    for fNo = 1:numel(flds)
        v = flds{fNo};
        if isfield(state2, v)    && ...
           isnumeric(state2.(v)) && ...
           all(size(state2.(v)) == size(state1.(v)))
           stateDiff.(v) = abs(state1.(v) - state2.(v));
           if opt.relative
               stateDiff.(v) = stateDiff.(v)./state1.(v);
           end
        else
            stateDiff = rmfield(stateDiff, v);
        end
    end
       
end
