function is_batch = isBatchMode()
% Check if MRST is running in batch mode
%
% SYNOPSIS:
%   is_batch = isBatchMode()
%
% NOTE:
%   This function simply checks if the global MRST_BATCH exists, and is
%   true if it exists.

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

    global MRST_BATCH
    if isempty(MRST_BATCH)
        MRST_BATCH = false;
    end
    assert(islogical(MRST_BATCH), 'MRST_BATCH global must be a true or false');
    is_batch = MRST_BATCH;
end
