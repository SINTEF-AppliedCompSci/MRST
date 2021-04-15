function grdecl = refineGrdecl(grdecl_in, dim)
% 
%   Refine an Eclipse grid (`GRDECL file`) with a specified factor in each of
%   the three logical grid directions.
% 
%   Currently, the function 
%     - Refines the grid by subdividing each cell according to input
%       parameter `dim`. 
%     
%     - Updates the following keywords: `ACTNUM`, `PERMX`, `PERMY`, `PERMZ`,
%       `PORO`, `MULTX`, `MULTY`, `MULTZ`, `EQLNUM`, `FLUXNUM` and `NTG`.
% 
%     - Updates the `FAULTS` keyword.  Multipliers are not changed, which is
%       only correct for refinements in the z-direction.
%    
%   Function does not handle flow-based upscaling keywords (e.g. `MULTX`,
%   `MULTY` and `MULTZ`) for which there is no natural automated refinement
%   process. 
% 
%
%   NOTE:
%     This function is not fully tested and has only been used on a limited
%     subset of models. While potentially useful, it should be used with care
%     and results should be carefully examined.
%
% SYNOPSIS:
%   function grdecl = refineGrdecl(grdecl_in, dim)
%
% PARAMETERS:
%   grdecl_in - Eclipse grid to refine (grdecl)
%   dim       - refinement factor in each logical direction (3-component vector)
%
% RETURNS:
%   grdecl - refined grid
%
% SEE ALSO:
% refineDeck

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
   
   [xyz, zcorn] = grdeclXYZ(grdecl_in);
   nx = dim(1); ny = dim(2); nz = dim(3);

   %% refining in z direction
   if(nz > 1)
      dz = zcorn(:, :, 2:2:end) - zcorn(:, :, 1:2:end);
      slice = zcorn(:, :, 1:2:end);
      zcorn = zeros(size(zcorn, 1), size(zcorn, 2), size(zcorn, 3) * nz);
      for i = 0:nz - 1
         zcorn(:, :, 2 * i + 1:2 * nz:end) = slice + i * dz / nz;
         zcorn(:, :, 2 * i + 2:2 * nz:end) = slice + (i + 1) * dz / nz;
      end
   end

   %% refining in x direction
   if(nx > 1)
      dx = xyz(:, 2:end, :) - xyz(:, 1:end - 1, :);
      slice = xyz(:, 1:end - 1, :);
      endslice = xyz(:, end, :);
      xyz = zeros(6, (size(xyz, 2) - 1) * nx + 1, size(xyz, 3));
      for i = 1:nx
         xyz(:, i:nx:end - 1, :) = slice + (i - 1) * dx / nx;
      end
      xyz(:, end, :) = endslice;
      
      dz = zcorn(2:2:end, :, :) - zcorn(1:2:end, :, :);
      slice = zcorn(1:2:end, :, :);
      zcorn = zeros(size(zcorn, 1) * nx, size(zcorn, 2), size(zcorn, 3));
      for i = 0:nx - 1
         zcorn(2 * i + 1:2 * nx:end, :, :) = slice + i * dz / nx;
         zcorn(2 * i + 2:2 * nx:end, :, :) = slice + (i + 1) * dz / nx;
      end
   end

   %% refining in y direction
   if(ny > 1)
      dy = xyz(:, :, 2:end) - xyz(:, :, 1:end - 1);
      slice = xyz(:, :, 1:end - 1);
      endslice = xyz(:,:, end);
      xyz = zeros(6, size(xyz, 2), (size(xyz, 3) - 1) * ny + 1);
      for i = 1:ny
         xyz(:, :, i:ny:end - 1) = slice + (i - 1) * dy / ny;
      end
      xyz(:, :, end) = endslice;
      
      dz = zcorn(:, 2:2:end, :) - zcorn(:, 1:2:end, :);
      slice = zcorn(:, 1:2:end, :);
      zcorn = zeros(size(zcorn, 1), size(zcorn, 2) * ny, size(zcorn, 3));
      for i = 0:ny - 1
         zcorn(:, 2 * i + 1:2 * ny:end, :) = slice + i * dz / ny;
         zcorn(:, 2 * i + 2:2 * ny:end, :) = slice + (i + 1) * dz / ny;
      end
   end

   grdecl.cartDims = grdecl_in.cartDims .* [nx, ny, nz];
   grdecl.COORD = xyz(:);
   grdecl.ZCORN = zcorn(:);

   %% refine faults, if any
   if isfield(grdecl_in, 'FAULTS')
      grdecl.FAULTS = refine_faults(grdecl_in.FAULTS, dim);
      if isfield(grdecl_in,'MULTFLT')
         grdecl.MULTFLT = grdecl_in.MULTFLT;
      end
   end
   
   %% expanding and copying associated cell-based fields
   cartDims_in = grdecl_in.cartDims;
   
   if ~isfield(grdecl_in, 'ACTNUM')
      % Interpret no ACTNUM as all cells active.  This is consistent with
      % (e.g.) 'processGRDECL'.
      grdecl_in.ACTNUM = ones([prod(cartDims_in), 1]);
   end
   
   % expanding and copying cell based fields
   fields = {'ACTNUM', 'PERMX', 'PERMY', 'PERMZ', 'PORO', 'MULTZ', 'MULTX', ...
             'MULTY', 'EQLNUM', 'FLUXNUM', 'NTG'};

   for f = fields
      if isfield(grdecl_in, f{:})
         A = reshape(grdecl_in.(f{:}), cartDims_in);
         grdecl.(f{:}) = A(ceil((1:nx * cartDims_in(1)) / nx), ...
                           ceil((1:ny * cartDims_in(2)) / ny), ...
                           ceil((1:nz * cartDims_in(3)) / nz));
         grdecl.(f{:}) = grdecl.(f{:})(:);
      end
   end
end  
 
% ----------------------------------------------------------------------------
function faults = refine_faults(faults_in, refs)
   
   sdim = {'X','Y','Z'};
   dd = size(faults_in);
   faults = faults_in;
   for i = 1:dd(1)
      for j = 0:2
         if ~(sdim{j+1}==faults_in{i,8})
            faults{i,1+2*j+1} = (faults_in{i,1+2*j+1}-1) * refs(j+1)+1;
            faults{i,1+2*j+2} = (faults_in{i,1+2*j+2}-1) * refs(j+1)+1+(refs(j+1)-1);
         else
            faults{i,1+2*j+1} = (faults_in{i,1+2*j+1}-1)*refs(j+1)+1+refs(j+1)-1;
            faults{i,1+2*j+2} = (faults_in{i,1+2*j+2}-1)*refs(j+1)+1+refs(j+1)-1;
         end
      end
   end
end