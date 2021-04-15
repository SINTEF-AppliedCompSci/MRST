function reginx = getRegMap(val, REGNUM, REGINX, varargin)
% Get region mapping from REGNUM array
%
% SYNOPSIS:
%   reginx = getRegMap(val, REGNUM, REGINX)
%
% REQUIRED PARAMETERS:
%   val      - Lookup values that are to be used to evaluate region tables.
%              Only used for dimension checking.
%
%   REGNUM   - Array of region indicators per cell in the simulation grid.
%
%   REGINX   - Override for output parameter. Only minimal checking is done
%              in this case.
%
% OPTIONAL PARAMETERS:
%   key - Parameter description
%
% RETURNS:
%   reginx - Cell array of length equal to number of regions, where each
%            entry corresponds to a list of cells belonging to that region.
% NOTE:
%   This is considered an internal function. The interface may change at
%   any time, or the function may disappear into the void without warning.
%
% SEE ALSO:
%   interpReg

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

   opt = struct('cellInx', []);
   opt = merge_options(opt, varargin{:});
   nt  = numel(REGINX);
   if numel(REGNUM) == 1
       reginx = cell(REGNUM,1);
       N = numel(double(val));
       reginx(REGNUM) = {1:N};
       return;
   end

   if isempty(val)
      % Allow for empty val
      if isempty(opt.cellInx)
         N = numel(REGNUM);
      else
         N = numel(opt.cellInx);
      end
   else
      N = numelValue(val);
   end

   if isempty(opt.cellInx)
      if nt == 1
         reginx = { ':' };
      else
         % Entire domain.
         if N ~= numel(REGNUM)
            % Do not use NUMEL in case of ADI.
            error('Region reference for input undefined');
         end
         reginx = REGINX;
      end
   elseif isempty(REGNUM)
      % If REGNUM is empty, there is only one region for all cells.
      reginx = {(1:numel(opt.cellInx))'};
   else
      % Reference to (small) subset of all cells
      cellInx = opt.cellInx;
      regnum  = REGNUM(cellInx);
      if numel(cellInx) > 1
         if N ~= numel(cellInx)
            % Do not use NUMEL in case of ADI.
            error('Number of cell indices must be same as input values');
         end
         reginx = arrayfun(@(x) find(x == regnum), 1 : nt, ...
                           'UniformOutput', false);
      elseif numel(cellInx) == 1
         % Allow single input (for exploring single cell functions).
         reginx         = cell(1, nt);
         reginx(regnum) = {(1:N)'};
      else
         error('Got empty cellInx input. This is not happening...');
      end
   end
end
