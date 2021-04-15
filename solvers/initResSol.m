function resSol = initResSol(G, p0, s0, varargin)
%Initialise incompressible reservoir solution data structure
%
% SYNOPSIS:
%   state = initResSol(G, p0)
%   state = initResSol(G, p0, s0)
%
% PARAMETERS:
%   G  - Grid data structure.
%
%   p0 - Initial reservoir pressure.  Scalar or a G.cells.num-by-1 vector.
%
%   s0 - Initial reservoir saturation.  A 1-by-(number of phases) vector
%        or a (G.cells.num)-by-(number of phases) array.
%        Default value: s0 = 0 (single phase).
%
% RETURNS:
%   state - Initialized reservoir solution structure having fields
%             - pressure -- One scalar pressure value for each cell in 'G'.
%             - flux     -- One Darcy flux value for each face in 'G'.
%             - s        -- Phase saturations for all phases in each cell.
%
% NOTE:
%   In the case of a `G.cells.num` by `3` array of fluid saturations,
%   state.s, the columns are generally interpreted as
%
%          1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapour
%
%   Single pressures (`p0`) and initial phase saturations (`s0`) are repeated
%   uniformly for all grid cells.
%
%   The initial Darcy flux is zero throughout the reservoir.
%
% SEE ALSO:
%   `initWellSol`, `incompTPFA`.

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

   if nargin <= 2, s0 = 1; end

   [nc, nf] = deal(G.cells.num, G.faces.num);
   if isfield(G, 'nnc') && isfield(G.nnc, 'cells')
       % Expand the number of interfaces with the number of non-neighboring
       % interfaces
       nf = nf + size(G.nnc.cells, 1);
   end

   if size(s0, 1) == 1
      s0 = repmat(s0(1,:), [nc, 1]);
   elseif size(s0, 1) ~= nc
      error(msgid('InitSat:Inconsistent'), ...
           ['Initial saturation must either be 1-by-np ', ...
            'or (G.cells.num (=%d))-by-np.'], nc);
   end

   if numel(p0) == 1
      p0 = repmat(p0, [nc, 1]);
   end

   resSol = struct('pressure', p0,             ...
                   'flux',     zeros([nf, 1]), ...
                   's',        s0);

   if nargin == 4
      if isa(varargin{1}, 'double') && ...
              size(varargin{1}, 2) == size(resSol.s, 2)
         z = varargin{1};
         if size(z, 1) == 1
            z = repmat(z, [nc, 1]);
         elseif size(z, 1) ~= nc
            error(msgid('Mass:WrongSize'), ...
                 ['Surface volume: Expected G.cells.num (=%d) ', ...
                  'rows. Got %d.'], nc, size(z, 1));
         end
         resSol.z = z;
      else
         warning(msgid('Mass:Inconsistent'), ...
                ['Initial mass (surface volume) must be ', ...
                 '(SIZE(s,2) (=%d))-column DOUBLE array (ignored).'], ...
                 size(resSol.s, 2));
      end
   end
end
