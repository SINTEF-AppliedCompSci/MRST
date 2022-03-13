function bc_cg = convertBC2Coarse(bc, G, CG, nsub, sub, coarse_f, sgn)
%Convert fine-scale boundary conditions to coarse scale.
%
% SYNOPSIS:
%   bc_cg = convertBC2Coarse(bc, G, CG, nsub, sub, coarse_f, sgn)
%
% DESCRIPTION:
%   Converts the fine boundary condition structure to a coarse boundary
%   structure for use in coarse transport simulations.
%
% REQUIRED PARAMETERS:
%   bc        - Fine grid boundary condition structures, as defined by
%               function addBC.
%
%   G, CG     - Grid, and coarse grid data structures respectively.
%
%   nsub, sub - Fine-scale subfaces constituting individual coarse grid
%               faces.  represented in a condensed storage format.
%               Typically computed using function subFaces.
%
%   sgn       - Sign of the fine-grid subface flux to be multiplied with
%               the corresponding flux, before summed along the coarse
%               interface to give the coarse net flux.
%
%   coarse_f  - The coarse face number for each fine-grid subface.
%
% RETURNS:
%   bc_cg     - Valid boundary condition structure, restricted to coarse
%               grid boundary.
%
% NOTE:
%   This function is presently restricted to flux boundary conditions.
%
% SEE ALSO:
%   `convertSource2Coarse`, `signOfFineFacesOnCoarseFaces`, `subFaces`.

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


   is_flux = strcmpi(bc.type, 'flux');
   is_pressure = strcmpi(bc.type, 'pressure');

   if all(is_flux),

      fface = bc.face(is_flux);

      % Find index of coarse face to which fine-faces belong
      f2c = sparse(sub, 1, coarse_f);

      % Coarse boundary faces
      cf             = false([CG.faces.num, 1]);
      cf(f2c(fface)) = true;

      % Restructuring of fine-grid fluxes and fine-grid saturations.
      flux         = zeros([G.faces.num, 1]);
      flux(fface)  = bc.value(is_flux);

      sat          = zeros([G.faces.num, size(bc.sat,2)]);
      sat(fface,:) = bc.sat(is_flux,:);

      % Accumulate flux and saturation to the coarse fields
      accum = sparse(coarse_f, 1:numel(sub), 1) * ...
              [bsxfun(@times, sgn .* flux(sub), sat(sub,:)), flux(sub)];

      coarse_sat = bsxfun(@rdivide, accum(:,1:end-1), accum(:,end));

      % Add faces, fluxes and saturations into a bc structure

      bc_cg = addBC([], find(cf), 'flux', accum(cf,end), ...
                    'sat', abs(coarse_sat(cf, :)));

   elseif all(is_pressure)
       ff_cf=nan(CG.parent.faces.num,1);
       cfacesno  = rldecode([1:CG.faces.num]',diff(CG.faces.connPos));
       ff_cf(CG.faces.fconn)=cfacesno;
       cbcf=ff_cf(bc.face);
       % % only use one value to deside boundary conditon
       [cbf,ci,cj]=unique(cbcf);
       bc_cg=struct('face',cbf,...
                    'type',{bc.type(ci)},...
                    'value',bc.value(ci));
   else
      error('convertBC2Coarse: only flux boundary conditions supported');
      bc_cg = []; return
   end
end
