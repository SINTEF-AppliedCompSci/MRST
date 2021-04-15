function mt_struct = addToTimeStruct(mt_struct, t_struct)
%Append single time step report structure to cumulative report structure
%
% SYNOPSIS:
%   allsteps = addToTimeStruct(allsteps, step)
%
% PARAMETERS:
%   allsteps - Report structure for all existing time steps.
%
%   step     - Report structure for a single time step.  Typically the
%              output of function 'wellCalculateProduction'.
%
% RETURNS:
%   allsteps - Report structure for all existing time steps, with 'step'
%              data appended (by means of HORZCAT).
%
% NOTE:
%   Any structure fields which exist in 'step' but not in the input
%   'allsteps' will be created in the output 'allsteps'.
%
% SEE ALSO:
%   `wellCalculateProduction`, `horzcat`.

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


   for f = reshape(fieldnames(t_struct), 1, []),
      nm = f{1};
      if ~isfield(mt_struct, nm), mt_struct.(nm) = []; end
      mt_struct.(nm) = [mt_struct.(nm), t_struct.(nm)];
   end
end
