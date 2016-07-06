function [p, z, fp, s] = impesAssembleStateVarsURC(state, bc, wells, wdp)
%Collect pressure and surface volumes into linear arrays
%
% SYNOPSIS:
%   [p, z, fp] = impesAssembleStateVars(state, bc, wells, wdp)
%
% PARAMETERS:
%   state - Reservoir and well solution structure.
%
%   bc    - Boundary condition data structure as defined by function
%           'addBC'.  If ISEMPTY(bc), then the face mobility on outer faces
%           will be that of the connecting cell.  Otherwise, the mobility
%           is taken from the boundary condition if it is an in-flow
%           condition and from the cell if it is an out-flow interface.
%
%   wells - Well data structure as defined by function 'addWell'.  If
%           ISEMPTY(wells), then no additional phase mobility will be
%           defined.  Otherwise the same provisions apply to 'wells' as to
%           the boundary condition structure 'bc'.
%
%   wdp   - Perforation gravity pressure adjustments.  This is the result
%           of modelling the behaviour of the free-stream flow within the
%           well track in the presence or absence of gravity effects.
%
%           Function 'impesTPFADefaultWellModel' is a simplified well model
%           that may be used to compute such values, but more advanced
%           approaches may be employed too.
%
% SEE ALSO:
%   impesComputeMobility, tpfaUpwindStateVars, impesTPFADefaultWellModel.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


   p  = state.pressure;
   z  = state.z;
   s  = state.s;
   fp = state.facePressure;

   if ~isempty(bc),
      di       = strcmpi(bc.type, 'pressure');
      bcp      = zeros([numel(bc.face), 1]);
      bcp( di) = bc.value(di);
      bcp(~di) = state.facePressure(bc.face(~di));

      p  = [p ; bcp   ];
      z  = [z ; bc.sat];
      s  = [s ; bc.sat];

      fp(bc.face(di)) = bcp(di);
   end

   if ~isempty(wells),
      nperf      = cellfun('prodofsize', { wells.cells });
      is_bhp     = strcmpi('bhp'       , { wells.type  });

      zw         = vertcat(wells.compi           );
      sw         = vertcat(wells.compi           );
      pw         = vertcat(state.wellSol.pressure);
      pw(is_bhp) = vertcat(wells(is_bhp).val     );

      zw = rldecode(zw, nperf);
      pw = rldecode(pw, nperf);

      assert (all(size(pw) == size(wdp)), ...
             ['Inconsistent definition of perforation gravity ', ...
              'pressure adjustment.']);

      p  = [p ; pw + wdp];
      z  = [z ; zw      ];
      s  = [s ; sw      ];
      fp = [fp; pw + wdp];
   end
end
