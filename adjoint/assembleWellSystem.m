function W = assembleWellSystem(G, W, varargin)
%Generate pressure linear system components for wells.
%
% SYNOPSIS:
%   W = assembleWellSystem(G, W)
%   W = assembleWellSystem(G, W, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G - Grid structure as described by grid_structure.
%
%   W - Well structure as defined by addWell &c.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Type -- The kind of system to assemble.  The choice made
%                         for this option influences which pressure solvers
%                         can be employed later.
%                         String.  Default value = 'hybrid'.
%                 Supported values are:
%                   - 'hybrid' : Hybrid system with inverse mass matrix
%                                (BI) for Schur complement reduction.
%                   - 'mixed'  : Hybrid system with regular mass matrix (B)
%                                for reduction to mixed form.
%
%                 Note: The value of this option should match the option
%                 pair passed to function computeMimeticIP.
%
% RETURNS:
%   W - Updated well structure.  Function assembleWellSystem adds field 'S'
%       to each well in 'W'.  The well system W(k).S has the following
%       fields:
%          S.BI/S.B -- Well hybrid system (inverse) 'mass' matrix.
%          S.C      -- Well hybrid system discrete divergence operator.
%          S.D      -- Well hybrid system flux continuity matrix.
%          S.RHS    -- Well hybrid linear system right hand side.
%
% SEE ALSO:
%   `computeMimeticIP`, `addWell`.

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


opt = struct('Type', 'hybrid');
opt = merge_options(opt, varargin{:});
systemType = opt.Type;

numWells = length(W);
for k = 1 : numWells,
   w = W(k);
   sizeB = [1, 1] * numel(w.cells);
   sizeC = [sizeB(1), G.cells.num];

   if ~strcmp(systemType, 'mixed'),
      S.BI = spdiags(     w.WI, 0, sizeB(1), sizeB(2));
   end
   if ~strcmp(systemType, 'hybrid')
      S.B  = spdiags(1 ./ w.WI, 0, sizeB(1), sizeB(2));
   end

   C = sparse(1:sizeC(1), w.cells, 1, sizeC(1), sizeC(2));
   D = sparse(ones([sizeB(1), 1]));
   if strcmp(w.type, 'bhp'),
      RHS.f = - w.val(ones([sizeB(1), 1]));
      RHS.h = 0;                % never used
   elseif strcmp(w.type, 'rate'),
      RHS.f = zeros([sizeB(1), 1]);
      RHS.h = - w.val;
   end

   S.sizeB = sizeB; S.sizeC = sizeC; S.sizeD = size(D);
   S.C = C; S.D = D;
   S.RHS = RHS;
   S.type = systemType;

   W(k).S = S;
end
