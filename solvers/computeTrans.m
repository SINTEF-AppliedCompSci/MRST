function T = computeTrans(G, rock, varargin)
%Compute transmissibilities.
%
% SYNOPSIS:
%   T = computeTrans(G, rock)
%   T = computeTrans(G, rock, 'pn', pv, ...)
%
% PARAMETERS:
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field `perm`.  The permeability
%          is assumed to be in measured in units of metres squared (m^2).
%          Use function `darcy` to convert from darcies to m^2, e.g::
%
%                 perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns for a
%          symmetric full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor::
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%
% OPTIONAL PARAMETERS:
%
%   K_system - Define the system permeability is defined in valid values
%              are `xyz` and `loc_xyz`.
%
%   cellCenters - Compute transmissibilities based on supplied
%                 `cellCenters` rather than default `G.cells.centroids`
%
%   cellFaceCenters - Compute transmissibilities based on supplied
%                     `cellFaceCenters` rather then default
%                     `G.faces.centroids(G.cells.faces(:,1), :)`
%
%   fixNegative     - Take absolute value of negative transmissibilities if
%                     present. This typically happens if the grid cell
%                     geometry is degenerate, for instance with a centroid
%                     outside the cell. Default: true.
%
% RETURNS:
%   T - half-transmissibilities for each local face of each grid cell in
%       the grid.  The number of half-transmissibilities equals the number
%       of rows in `G.cells.faces`.
%
% NOTE:
%   Face normals are assumed to have length equal to the corresponding face
%   areas. This property is guaranteed by function `computeGeometry`.
%
% SEE ALSO:
%   `computeGeometry`, `permTensor`, `makeRock`.

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

   opt = struct('verbose'        , false, ...
                'grdecl'         , []   , ...
                'K_system'       , 'xyz', ...
                'cellCenters'    , []   , ...
                'fixNegative'    , true, ...
                'cellFaceCenters', []   );
   opt = merge_options(opt, varargin{:});

   dispif(opt.verbose, 'Computing one-sided transmissibilities...\t');
   t0 = ticif(opt.verbose);

   switch lower(opt.K_system)
      case 'xyz',     T = htrans_xyz  (G, rock, opt);
      case 'loc_xyz', T = htrans_local(G, rock, opt);

      otherwise
         error(msgid('PermSys:Unknown'), ...
               'Unknown permeability coordinate system ''%s''.', ...
               opt.K_system);
   end
   if opt.fixNegative
      T = handle_negative_trans(T, opt);
   end
   T = handle_ntg(T, rock, G);
   
   if isstruct(opt.grdecl) && numel(opt.grdecl) == 1
      m = getMultipliers(G, opt.grdecl);
      T = T .* m;
   end

   tocif(opt.verbose, t0);
end

%--------------------------------------------------------------------------

function T = htrans_xyz(G, rock, opt)
   [C, N, cellNo] = cell_geometry(G, opt);
   [K, i, j] = permTensor(rock, G.griddim);

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);

   % Compute T = C'*K*N / C'*C. Loop-based to limit memory use.
   T = zeros(size(cellNo));
   for k = 1 : size(i, 2)
      T = T + (C(:, i(k)) .* K(cellNo, k) .* N(:, j(k)));
   end

   T = T ./ sum(C .* C, 2);
end

%--------------------------------------------------------------------------

function T = htrans_local(G, rock, opt)
   if size(rock.perm, 2) == 1
      rock.perm = repmat(rock.perm, [1, G.griddim]);
   end

   if size(rock.perm, 2) ~= size(G.cartDims, 2)
      error(msgid('PermSys:Inconsistent'),                 ...
           ['Permeability coordinate system ''loc_xyz'' ', ...
            'is only valid for diagonal tensors.']);
   end

   assert (size(rock.perm, 1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(rock.perm,1), G.cells.num);

   [C, N, cellNo] = cell_geometry(G, opt);

   dim = ceil(G.cells.faces(:,2) / 2);
   ind = sub2ind(size(rock.perm), cellNo, double(dim));
   T   = reshape(rock.perm(ind), [], 1) .* ...
         sum(C .* N, 2) ./ sum(C .* C, 2);
     
end

%--------------------------------------------------------------------------

function [C, N, cellNo] = cell_geometry(G, opt)
   % Vectors from cell centroids to face centroids
   cellNo = gridCellNo(G);

   if ~isempty(opt.cellCenters)
      C = opt.cellCenters;
   else
      C = G.cells.centroids;
   end

   if ~isempty(opt.cellFaceCenters)
      C = opt.cellFaceCenters - C(cellNo,:);
   else
      C = G.faces.centroids(G.cells.faces(:,1), :) - C(cellNo,:);
   end

   % Outward-pointing normal vectors
   cf  = G.cells.faces(:,1);
   sgn = 2*(cellNo == G.faces.neighbors(cf, 1)) - 1;
   N   = bsxfun(@times, sgn, G.faces.normals(cf, :));
end

%--------------------------------------------------------------------------

function T = handle_ntg(T, rock, G)
if isfield(rock, 'ntg') && size(G.cells.faces, 2) == 2
    cellNo = gridCellNo(G);
    dim    = ceil(G.cells.faces(:,2) / 2);
    m      = rock.ntg(cellNo);
    m(dim==3) = 1;
    T = m.*T;
end
end
    
%--------------------------------------------------------------------------

function T = handle_negative_trans(T, opt)
   is_neg = T < 0;

   if any(is_neg)
      dispif(opt.verbose, ...
            ['\nWarning:\n\t%d negative transmissibilities.\n\t', ...
             'Replaced by absolute values...\n'], sum(is_neg));

      T(is_neg) = -T(is_neg);
   end
end

%--------------------------------------------------------------------------

function multipliers = getMultipliers(G, grdecl)
%Extract multipliers from grdecl.FAULTS and grdecl.MULT*
%
% SYNOPSIS:
%   m = getMultipliers(G, rock)
%
% PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   grdecl  - Raw ECLIPSE data field as read by (???).  The fields FAULTS
%             MULTX, MULTY and MULTZ all contain transmissibility
%             multipliers.  These are interpreted as one-sided multipliers,
%             i.e multiplicative modifiers for the half-face
%             transmissibilities computed in this file.
%
% RETURNS:
%   m - half-transmissibilities multipliers for all faces.

   multipliers = ones(size(G.cells.faces, 1), 1);

   % do fault multiplicators which are in grdecl.FAULTS.('fault_name')
   dirstruct = struct('X', 1, 'Y', 3, 'Z', 5);
   if isfield(grdecl,'FAULTS')
      faults = struct2cell(grdecl.FAULTS);
      for i = 1:numel(faults)

         fault = faults{i};

         % Skip if there are no fault multipliers.
         if ~isfield(fault, 'multflt')
            continue;
         end

         for j = 1 : numel(fault.dir)
            dirnum = dirstruct.(fault.dir(j));
            region = fault.cells(j,:);

            if region(dirnum) ~= region(dirnum+1)
               error(['Wrong fault ind direction ', fault.dir(j)]);
            end

            % For each half-face...
            for side = 0 : 1
               % Find cell indices in G corresponding to region
               region(dirnum:dirnum+1) = region(dirnum:dirnum+1) + side;

               [X,Y,Z] = ndgrid(region(1) : region(2), ...
                                region(3) : region(4), ...
                                region(5) : region(6));
               ind = sub2ind(G.cartDims, X(:), Y(:), Z(:));
               c   = cart2active(G, ind);

               % Set multiplier of cellFace with the correct tag
               cf  = mcolon(G.cells.facePos(c), G.cells.facePos(c+1)-1);
               tag = dirnum + 1 - side;
               cf  = cf(G.cells.faces(cf,2) == tag);
               multipliers(cf) = fault.multflt;
            end
         end
      end
   end

   %assume MULT?M is not used
   vecstruct=[];
   dirstruct=[];
   if isfield(grdecl,'MULTX')
      vecstruct = [vecstruct, grdecl.MULTX];
      dirstruct = [dirstruct, 2];
      mult      = ones(G.cartDims);
      multx     = reshape(grdecl.MULTX, G.cartDims);
      mult(2:end,:,:) = multx(1:end-1, :, :);
      vecstruct = [vecstruct, mult(:)];
      dirstruct = [dirstruct, 1];%set the direction number
   end

   if isfield(grdecl,'MULTY')
      vecstruct = [vecstruct, grdecl.MULTY];
      dirstruct = [dirstruct, 4];
      mult      = ones(G.cartDims);
      multx     = reshape(grdecl.MULTY, G.cartDims);
      mult(:,2:end,:) = multx(:, 1:end-1, :);
      vecstruct = [vecstruct, mult(:)];
      dirstruct = [dirstruct, 3];
   end

   if isfield(grdecl,'MULTZ')
      vecstruct = [vecstruct, grdecl.MULTZ];
      dirstruct = [dirstruct, 6];
      mult      = ones(G.cartDims);
      multx     = reshape(grdecl.MULTZ, G.cartDims);
      mult(:,:,2:end) = multx(:, :, 1:end-1);
      vecstruct = [vecstruct, mult(:)];
      dirstruct = [dirstruct, 5];
   end

   if ~isempty(vecstruct)
      % do multipliers
      for i = 1:length(dirstruct)
         vec   = vecstruct(:, i);
         dir   = dirstruct(i);
         b     = cumsum(G.cells.faces(:, 2) == dir);
         cfind = cumsum(double(diff(G.cells.facePos)));%facepos
         numf  = b(cfind)-[0; b(cfind(1:end-1))];
         ci    = G.cells.faces(:,2) == dirstruct(i);
         mult  = rldecode(vec(G.cells.indexMap), numf);

         multipliers(ci) = mult;
      end
   end
end
