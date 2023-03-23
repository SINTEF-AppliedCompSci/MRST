function [T,hybrids] = computeTrans_DFM(G, rock, varargin)
% Compute transmissibilities using a two-point scheme.
%
% SYNOPSIS:
%   T = computeTrans_DFM(G, rock)
%   T = computeTrans_DFM(G, rock, 'pn', pv, ...)
%
% This function is modified from the original computeTrans to include
% hybrid (essentially lower-dimensional) cells representing fractures etc.
% The hybrid cells are located on the interface between two 'normal' cells
% and have no extension (they are lines in 2D and planes in 3D) in the
% geometric grid.
% However, they are accounted for in the computational grid and are
% assigned an area (2D) or volume (3D). For more information see
% Karimi-Fard (SPE 2004).
%
% This function is a part of the DFM-module of MRST, and should be applied
% together with other functions of that module.
%
% PARAMETERS:
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field 'perm'.  The permeability
%          is assumed to be in measured in units of meters squared (m^2).
%          Use function 'darcy' to convert from darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns for a
%          symmetric full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor
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
%   K_system - Define the system permeability is defined in valid values
%              are 'xyz' and 'loc_xyz'.
%   cellCenters - Compute transmissibilities based on supplied cellCenters
%              rather than default G.cells.centroids
%   cellFaceCenters - Compute transmissibilities based on supplied
%              cellFaceCenters rather then default
%              G.faces.centroids(G.cells.faces(:,1), :)
%   hybrid - if true the distance between the cell center and the face center
%              is extended half an aperture in the direction of the normal.
%   compactNeighboringCell - Whether hybrid faces should be expanded for
%                            the cells on both sides of the face. Only
%                            activated if hybrid is true
%
% RETURNS:
%   T - half-transmissibilities for each local face of each grid cell in
%       the grid.  The number of half-transmissibilities equals the number
%       of rows in G.cells.faces.
%
% COMMENTS:
%   PLEASE NOTE: Face normals are assumed to have length equal to the
%   corresponding face areas.  This property is guaranteed by function
%   'computeGeometry'.
%
% SEE ALSO:
%   `computeGeometry`, `computeMimeticIP`, `darcy`, `permTensor`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen. 2013 Iris AS

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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


   opt = struct('verbose', false, 'grdecl', [], 'K_system', 'xyz', ...
                'cellCenters', [], 'cellFaceCenters', [],'hybrid',false, ...
                'compactNeighboringCell',false);
   opt = merge_options(opt, varargin{:});

   if ~ischar(opt.K_system),
      error(msgid('PermSys:UnsupportedType'), ...
           ['Specified permeability coordinate system must be a ', ...
            'string value.']);
   end

   dispif(opt.verbose, 'Computing one-sided transmissibilities...\t');
   t0 = ticif(opt.verbose);

   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   sgn    = 2*(cellNo == G.faces.neighbors(G.cells.faces(:,1), 1)) - 1;
   N      = bsxfun(@times, sgn, G.faces.normals(G.cells.faces(:,1),:));

   if ~isempty(opt.cellCenters)
      cellCenters = opt.cellCenters;
   else
      cellCenters = G.cells.centroids;
   end

   if ~isempty(opt.cellFaceCenters)
      cellFaceCenters = opt.cellFaceCenters;
   else
      cellFaceCenters = G.faces.centroids(G.cells.faces(:,1), :);
   end

   C      = cellFaceCenters - cellCenters(cellNo, :);

   % Faces between hybrid and non-hybrid cells must have a centroid that is
   % different from the center of the hybrid cell.
   if opt.hybrid

       % Find faces on the interface between hybrid and non-hybrid cells
       hybrids = (G.cells.hybrid(cellNo)==1 & ~G.faces.hybrid(G.cells.faces));

       % There should be more than one such face (unless the is only one
       % hybrid cell, which moreover is located on the boundary of the
       % domain).
       if sum(hybrids) > 1

           % Move the location of the face centroid (for this computation
           % only) with half the aperture of the hybrid cell in the
           % direction of the normal vector of the face
           n = N(hybrids,:);

           % Aperture is volume / area of the face
           apertures = G.cells.volumes(cellNo(hybrids))./G.faces.areas(G.cells.faces(hybrids));

           % Account for the area-weighting of the normal vectors
           expantion = bsxfun(@times,n,apertures./2./G.faces.areas(G.cells.faces(hybrids)));
           C(hybrids,:) = expantion;

           % By default there is a slight inconsistency in how the
           % centroids are defined; it is expanded as seen from the hybrid
           % cell, but keeps its coordinates as seen from the normal cell
           % on the other side. The latter is modified by the code
           % below.
           % Note that if the neighboring cell is small relative to the
           % aperture, this may cause a negative transmissibility.
           if opt.compactNeighboringCell
               % Find the index of the faces as seen from the 'normal'
               % cells on the other side of the connections
               next2hybrids = getNext2hybridInterfaces(G,hybrids);

               % Do the extension there as well
               C(next2hybrids,:) = C(next2hybrids,:) + expantion;
           end
       end

	   % Do a similar exantion for the hybrid cells of second kind
       if any(G.cells.hybrid == 2)
           hybrids2 = G.cells.hybrid(cellNo)==2;
           n = N(hybrids2,:);
           apertures = sqrt(G.cells.volumes(cellNo(hybrids2)));
           expantion = bsxfun(@times,n,apertures./2./G.faces.areas(G.cells.faces(hybrids2)));
           C(hybrids2,:) = expantion;
       end

   end

   C2     = sum(C.*C, 2);
   if strcmpi(opt.K_system, 'xyz'),
      [K, i, j] = permTensor(rock, size(G.nodes.coords,2));
      %T = sum(C(:,i) .* K(cellNo,:) .* C(:,j), 2) ./ C2 .* ...
      T = sum(C(:,i) .* K(cellNo,:) .* N(:,j), 2) ./ C2;
   elseif strcmpi(opt.K_system, 'loc_xyz'),

      if size(rock.perm,2) ~= size(G.cartDims,2),
         error(msgid('PermSys:Inconsistent'),                 ...
              ['Permeability coordinate system ''loc_xyz'' ', ...
               'is only valid for diagonal tensor']);
      end

      dim = ceil(G.cells.faces(:,2) / 2);
      ind = sub2ind(size(rock.perm), cellNo, double(dim));
      T   = reshape(rock.perm(ind), [], 1) .* sum(C .* N, 2) ./ C2;
   else
      error(msgid('PermSys:Unknown'), ...
            'Unknown permeability coordinate system ''%s''.', ...
            opt.K_system);
   end

   is_neg = T < 0;
   if any(is_neg),
      dispif(opt.verbose, ...
            ['\nWarning:\n\t%d negative transmissibilities.\n\t', ...
             'Replaced by absolute values...\n'], sum(is_neg));

      T(is_neg) = -T(is_neg);
   end
   clear is_neg

   if isstruct(opt.grdecl) && numel(opt.grdecl) == 1,
      m = getMultipliers(G, opt.grdecl);
      T = T .* m;
   end

   tocif(opt.verbose, t0);
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
   if isfield(grdecl,'FAULTS'),
      faults = struct2cell(grdecl.FAULTS);
      for i = 1:numel(faults),

         fault = faults{i};

         % Skip if there are no fault multipliers.
         if ~isfield(fault, 'multflt'),
            continue;
         end

         for j = 1 : numel(fault.dir),
            dirnum = dirstruct.(fault.dir(j));
            region = fault.cells(j,:);

            if region(dirnum) ~= region(dirnum+1),
               error(['Wrong fault ind direction ', fault.dir(j)]);
            end

            % For each half-face...
            for side = 0 : 1,
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
   if isfield(grdecl,'MULTX'),
      vecstruct = [vecstruct, grdecl.MULTX];
      dirstruct = [dirstruct, 2];
      mult      = ones(G.cartDims);
      multx     = reshape(grdecl.MULTX, G.cartDims);
      mult(2:end,:,:) = multx(1:end-1, :, :);
      vecstruct = [vecstruct, mult(:)];
      dirstruct = [dirstruct, 1];%set the direction number
   end

   if isfield(grdecl,'MULTY'),
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
      for i = 1:length(dirstruct),
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

function there = getNext2hybridInterfaces(G,here)
% For a set of indexes (here) in terms of G.cells.faces, find the index of
% the face (in terms of G.cells.faces) from the neighboring cell.
% NOTE: The method is tailored for hybrid cells; it is tacitly assumed that
% here are between hybrid and normal cells

% Cell faces
fno = G.cells.faces;

% Find map between to the other cells.
% We now that the hybrid cells are stored at the end of G.cells.faces, so
% the occurence we are looking for is the first one.
[~,hereMap]=unique(fno,'first');

% pick out the other face.
there = hereMap(fno(here));

end
