function bc = coarsenBC(cg, bc)
%Construct coarse-grid boundary conditions from fine-grid boundary cond.
%
% SYNOPSIS:
%   bc = coarsenBC(cg, bc)
%
% DESCRIPTION:
%   All 'pressure' boundary conditions are mapped from fine-grid to the
%   coarse-grid faces by sampling the value at cg.faces.centerFace *). All
%   flux boundaries are accumulated
%
% PARAMETERS:
%   cg      - Coarse-grid including parent grid.
%
%   bc      - MRST boundary conditions struct as returned by addBC.
%
% RETURNS:
%   bc      - MRST boundary conditions for coarse-grid cg.  For 'pressure'
%             conditions, the value is sampled in cg.faces.centerFace.  For
%             'flux' conditions, the nex flux accoss the coarse face is
%             computed.  Mixed conditions are not supported.
%
% SEE ALSO:
%   `coarsenGeometry`, `coarsenFlux`, `fineToCoarseSign`, `upscaleSchedule`

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

   if isempty(bc),
         return;
   end

   assert(isfield(cg, 'parent'), ...
      'Huh!? Field ''parent''missing in coarse grid.  Did you really supply a coarse grid?');

   % Construct map from fine-to-coarse map for faces
   cf  = zeros(cg.parent.faces.num, 1);

   % Extract one face per coarse block to sample fine grid. This could be
   % improved by finding the actual center face instead of just the first.
   centerFace = cg.faces.fconn(cg.faces.connPos(1:end-1));
   cf(centerFace) = 1:cg.faces.num;

   % For 'flux' boundary conditions, account for possible sign change
   sgn   = nan(cg.parent.faces.num, 1);
   sgn(cg.faces.fconn)   = fineToCoarseSign(cg);

   ix    = cf > 0;
   face  = bc.face( ix(bc.face));
   v     = bc.value(ix(bc.face));
   type  = bc.type( ix(bc.face));

   ix    = strcmpi(type, 'flux');
   v(ix) = sgn(face(ix)) .* v(ix);

   %%% >>> HACK to accumulate flux boundary conditions (lacking sign change)
   fix=rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';

   ispressure = false(cg.parent.faces.num, 1);
   ispressure (bc.face) = strcmpi(bc.type, 'pressure');

   i = accumarray(fix, ispressure(cg.faces.fconn));
   j = (i == diff(cg.faces.connPos)) | i == 0;

   if ~all(j),
      error('Mixed boundry conditions on the coarse scale is unsupported.');
   end

   val = zeros(cg.parent.faces.num, 1);
   val(bc.face)=bc.value;
   val = accumarray(fix, val(cg.faces.fconn));
   val(i>0)=0;

   k = strcmpi(type, 'flux');

   v(k) = val(cf(face(k)));
   %%% <<<

   % Sample bc in fine faces present in the fine-to-coarse map "cf"
   bc = addBC([], cf(face), type, v);
end
