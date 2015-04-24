function varargout = impesComputeMobility(state, fluid, bc, wells, wdp)
%Compute phase mobility at dynamic state.
%
% SYNOPSIS:
%    mob        = impesComputeMobility(state, fluid, bc, wells, wdp)
%   [mob, dmob] = impesComputeMobility(...)
%
% PARAMETERS:
%   state - Reservoir and well solution structure.
%
%   fluid - Fluid object.
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
% RETURNS:
%   mob  - Phase mobility.  One scalar, non-negative value for each fluid
%          phase (at reservoir conditions) for each cell represented by
%          'state' and each perforation and each boundary condition if
%          applicable.
%
%   dmob - Phase mobility derivatives.  Specifically, the derivatives of
%          'mob' with respect to each phase saturation.  Ordered with
%          mobility phase index cycling the most rapidly, i.e.
%
%                           d mob(i,MOD(j-1,np)+1)
%               dmob(i,j) = ----------------------
%                              ds(i,CEIL(j/np))
%
%          in which 'np' denotes the number of fluid phases at reservoir
%          conditions.  In other words, the mobility Jacobian is stored in
%          column-major (Fortran) ordering.  We furthermore assume that the
%          number of fluid phases at reservoir conditions equals the number
%          of fluid components at surface conditions.
%
% NOTE:
%   The return values are ordered according to convention in the
%   'impesTPFA' solver and related functions.
%
% SEE ALSO:
%   impesAssembleStateVars, tpfaUpwindStateVars, impesTPFADefaultWellModel.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   [p, z, fp, s] = impesAssembleStateVarsURC(state, bc, wells, wdp);

   if ~all(p > 0),
      error('Unphysical non-positive pressure values detected.');
   end

   if any(any(z < 0)),
      error('Unphysical negative mass values detected.');
   end

   [mu, mu, mu, u] = fluid.pvt(p, z);                                  %#ok

   if ~all(all(mu > 0)),
      error('Unphysical non-positive viscosity values detected.');
   end

   if any(any(u < -sqrt(eps))),
      error('Unphysical negative phase volume values detected.');
   end

   %s               = state.s;

   [kr{1:nargout}] = fluid.relperm(s);

   if any(any(kr{1} < -sqrt(eps))),
      error('Unphysical negative rel-perm values detected.');
   end

   varargout{1}    = kr{1} ./ mu;

   if nargout > 1,
      % Caller requested mobility derivatives too...
      %
      varargout{2} = kr{2} ./ repmat(mu, [1, size(kr{1}, 2)]);
   end
end
