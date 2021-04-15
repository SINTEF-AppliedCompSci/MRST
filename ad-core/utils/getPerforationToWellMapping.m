function [perf2well, Rw] = getPerforationToWellMapping(w)
%Get map from global perforation number to global well index.
%
% SYNOPSIS:
%   perf2well = getPerforationToWellMapping(w)
%
% REQUIRED PARAMETERS:
%   w          - Well structure.
%
% RETURNS:
%   perf2well  - perf2well(ix) will give the well number of global
%                perforation number ix.
%
% SEE ALSO:
%   WellModel

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
    if isempty(w)
        perf2well = [];
        Rw = [];
        return
    end
    nw = numel(w);
    nConn       = cellfun(@numel, {w.cells})'; % # connections of each well
    perf2well   = rldecode((1:nw)', nConn);
    if nargout > 1
        nperf = numel(perf2well);
        if nperf == nw
            Rw = 1;
        else
            Rw = sparse((1:nperf)', perf2well, 1, nperf, nw);
        end
    end
end
