function state = initState(G, W, p0, varargin)
%Initialise state object for reservoir and wells.
%
% SYNOPSIS:
%   state = initState(G, W, p0)
%   state = initState(G, W, p0, s0)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   W  - Well structure.  Pass an empty array for models without wells.
%
%   p0 - Initial reservoir pressure.  Also used as initial well BHP.
%
%   s0 - Initial reservoir composition (saturation).  If the model contains
%        wells, then the number of declared components in the injected well
%        fluids must equal the number of declared components in the
%        reservoir fluids.
%
% RETURNS:
%   state - Initial reservoir (and well) state object.
%
% SEE ALSO:
%   `grid_structure`, `addWell`.

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


   state = initResSol(G, p0, varargin{:});

   if ~isempty(W)
      state.wellSol = initWellSol(W, p0);

      if nargin > 3
         s0 = varargin{1};

         try
            compi = vertcat(W.compi);
         catch
            error('Well compositions are inconsistently specified');
         end

         assert (size(compi, 2) == size(s0, 2),                        ...
                ['Number of phases in well definition (%d) does not ', ...
                 'match number of phases in reservoir fluids (%d)\n' , ...
                 'Consider using the ''Comp_i'' option of '          , ...
                 'function ''addWell''.'], size(compi, 2), size(s0, 2));
      end
   end
end
