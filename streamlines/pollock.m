function varargout = pollock(G, state, varargin)
%Trace streamlines in logically Cartesian grid using Pollock approximation.
%
% SYNOPSIS:
%   [S,T,C] = pollock(G, state)
%   [S,T,C] = pollock(G, state, startpos)
%   [S,T,C] = pollock(G, state, 'pn', pv, ...)
%   [S,T,C] = pollock(G, state, startpos, 'pn', pv, ...)
%
% PARAMETERS:
%
%   G         - Cartesian or logically Cartesian grid.
%
%   state     - State structure with field 'flux'.
%
% OPTIONAL PARAMETER
%
%   positions - Matrix of size (N, 1) or (N, d+1), where d is the dimension
%               of the grid, used to indicate where the streamlines should
%               start.
%
%               If the size is (N, 1), positions contains the cell indices
%               in which streamlines should start. Each streamline is
%               started in the the local coordinate (0.5, 0.5, ...). To be
%               precise, this is the mean of the corner points, not the
%               centroid of the cell.
%
%               If the size is (N, d+1), the first column contains cell
%               indices, and the d next columns contain the local
%               coordinates at which to start streamlines.
%
% OPTIONAL PARAMETERS:
%
%   substeps  - Number of substeps in each cell, to improve visual quality.
%               Default 5.
%
%   maxsteps  - Maximal number of points in a streamline.
%               Default 1000.
%
%   reverse   - Reverse velocity field before tracing.
%               Default false.
%
%   pvol      - Pore-volume vector.  One positive scalar for each active
%               cell in the grid `G`.  Makes the physical interpretation of
%               time-of-flight appropriate for non-uniform porosity fields.
%
%   isoutflow - Cell-wise boolean indicator which indicates if a streamline
%               should terminate upon reaching that cell. Defaults to
%               false(G.cells.num, 1) and is useful in the presence of many
%               weak source terms (which do not lead to inflow or outflow
%               over all faces for a given source term).
%
%   blocksize - Internal parameter indicating how many streamlines are
%               processed simultaneously. Larger values give faster
%               processing, at a higher memory cost. Default: 1000.
%
% RETURNS:
%
%  S      - Cell array of individual streamlines suitable for calls like
%           streamline(pollock(...)) and streamtube(pollock(...)).
%
%  T      - Time-of-flight of coordinate.
%
%  C      - Cell number of streamline segment, i.e, line segment between
%           two streamline coordinates.
%
% EXAMPLE:
%
%   S = pollock(G, x);
%   % pad with nan's
%   S = reshape([S, repmat({[nan, nan]}, [numel(S),1])]',[], 1);
%   S = vertcat(S{:});
%   plot(S(:,1), S(:,2), 'r-');
%
%   streamline(pollock(G, x, 'pvol', poreVolume(G, rock)));
%
% SEE ALSO:
%   `streamline`.

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

% Written by Jostein R. Natvig, SINTEF ICT, 2010.

   d = size(G.nodes.coords, 2);
   if mod(nargin, 2) == 0
      positions = [(1:G.cells.num)', repmat(0.5, [G.cells.num, d])];
   else
      positions = varargin{1};
      if size(positions, 2) == 1
         positions = [positions, repmat(0.5, [size(positions, 1), d])];
      elseif size(positions, 2) ~= 1 + d
         error('Expected array of local positions of width 1 or 1+d.');
      end
      varargin = varargin(2:end);
   end

   opt = struct('substeps'  , 5,     ...
                'maxsteps'  , 1000,  ...
                'reverse'   , false, ...
                'isoutflow' , false(G.cells.num, 1), ...
                'blocksize' , 1000,  ...
                'pvol'      , ones([G.cells.num, 1]), ...
                'flux'      , [], ... % See getInteriorFluxesStreamlines
                'neighbors' , []);    % See findStreamlineNeighborshipCart
   opt = merge_options(opt, varargin{:});

   if ~all(opt.pvol > 0)
      error('PoreVol:NonPositive', ...
            'Rock pore-volume must be strictly positive in all cells');
   end

   [varargout{1:nargout}] = trace(G, state, positions, opt);
end

% =========================================================================

function varargout = trace(G, state, pos, opt)
   d              = size(G.nodes.coords, 2);
   numStreamlines = size(pos, 1);
   assert(size(pos, 2) == d + 1);

   if ~isfield(G, 'cellNodes')
      cn = cellNodes(G);
      G.cellNodes = accumarray(cn(:,1:2), cn(:,3));
   end

   % Make array face fluxes for each cell in grid (Not outer).
   flux = opt.flux;
   if isempty(flux)
      flux = getInteriorFluxesStreamlines(G, state, opt.pvol, opt.reverse);
   end

   neighbors = opt.neighbors;
   if isempty(neighbors)
      neighbors  = findStreamlineNeighborshipCart(G);
   end

   magic  = opt.blocksize;
   XYZ    = nan(numStreamlines, d, magic);
   T      = nan(numStreamlines, magic);
   C      = nan(numStreamlines, magic);
   active = true(numStreamlines, 1);


   % Store crossing coordinates of active streamlines
   [XYZ(active,:,1)] = globalCoordinate(G, pos(active,1), pos(active, 2:end));
   T(active, 1) = zeros(sum(active), 1);
   C(active, 1) = pos(active,1);

   i = 2;
   while any(active)
      % Realloc
      if i+opt.substeps+1 > size(XYZ, 3)
         magic = max(magic, opt.substeps+1);
         XYZ   = cat(3, XYZ, nan(numStreamlines, d, magic));
         T     = cat(2, T,   nan(numStreamlines, magic));
         C     = cat(2, C,   nan(numStreamlines, magic));
      end
      current_cell = pos(active,1);

      % Take another pollock step
      [pos(active, :), t, xyz] = step(pos(active,:), flux, neighbors, opt.substeps);

      % Store crossing coordinates and, optionally, coordinates along curve
      % trajectory in cell of active streamlines
      for k=1:opt.substeps
         [XYZ(active, :, i+k-1)] = globalCoordinate(G, current_cell, xyz(:,:,k));
      end
      T(active, i-1+(1:opt.substeps)) = repmat(t/opt.substeps, [1, opt.substeps]);
      C(active, i-1+(1:opt.substeps)) = repmat(pos(active, 1), [1, opt.substeps]);

      % Update active flag
      active(active)    =  pos(active,1) ~= current_cell & ~opt.isoutflow(current_cell);

      i = i+opt.substeps;
      if i > opt.maxsteps
          break;
      end
   end

   % Pack coordinates in list with streamlines separated by NaN.
   p = reshape(permute(XYZ, [3,1,2]), [], d);

   i = ~isnan(p(:,1));
   j = i|[true;i(1:end-1)];
   p = p(j,:);

   % Pack streamline coordinates in a cell array suitable for use with
   % Matlab streamline, i.e., as in 'streamline(pollock(G, resSol));'
   flag = isnan(p(:,1));
   ix = find(flag);
   dd  = diff([0;ix])-1;
   varargout{1} = mat2cell(p(~flag,:), dd, d);
   if nargout > 1
      T = reshape(T', [], 1);
      T = T(j);
      varargout{2} = mat2cell(T(~flag), dd, 1);
   end
   if nargout > 2
      C = reshape(C', [], 1);
      C = C(j);
      varargout{3} = mat2cell(C(~flag), dd, 1);
   end
end

% =========================================================================

function xyz = globalCoordinate(G, c, p)
% Compute global coordinate corresponding to local coordinate p in cells c
% p  - local positions == [xi,eta,zeta] in 3D
% c  -

   if numel(c)==1, p = reshape(p, 1, []); end

   % Compute node weight for quadrilateral or hexahedron
   d = size(G.nodes.coords, 2);
   w = ones(size(p,1), 2^d);
   for i=1:d
      mask        = logical(bitget((0:2^d-1)', i));
      w(:, mask)  = w(:, mask).* repmat(  p(:,i), [1, sum( mask)]);
      w(:,~mask)  = w(:,~mask).* repmat(1-p(:,i), [1, sum(~mask)]);
   end

   % Compute weighted average of corner points
   xyz = zeros(size(p,1), d);
   for i=1:d
      xi       = G.nodes.coords(:,i);
      xyz(:,i) = sum(w .* reshape(xi(G.cellNodes(c, :))', 2^d, [])', 2);
   end
end

% =========================================================================

function [pos, tof, xyz] = step(pos, flux, neighbors, nsubsteps)
% Update pos array by computing new local coordinate and new cell.
% In addition, compute curve within cell.

   f = flux(pos(:,1),:);
   n = neighbors(pos(:,1),:);

   dims = size(pos, 2)-1;
   T    = nan(size(pos,1),dims);
   for i=1:dims
      T(:,i) = computeTime(pos(:,1+i), f(:,2*i-1:2*i));
   end
   [tof, dir] = min(T, [], 2);

   xyz = zeros(size(pos,1), dims, nsubsteps);
   d   = zeros(size(pos, 1), 1);
   for s=1:nsubsteps
      for i=1:dims
         t = tof*s/nsubsteps;
         [xyz(:,i,s), d(:,i)] = computePosition(pos(:,1+i), f(:,2*i-1:2*i), t);
      end
   end

   pos (:,2:end) = xyz(:,:,s);

   % Find direction to look up neighbor cell
   k  = 2*(dir-1)+d(sub2ind([numel(dir), 3], (1:numel(dir))', dir));
   t  = sub2ind(size(n), (1:numel(k))', k);

   % Update cell number if NOT at boundary.
   % IF at boundary, mark dir with NaN to avoid changing local coordinate
   % below.
   ind         = n(t)==0;
   % Also, if there is no finite time to escape the current cell, set NaN
   % value to indicate that the streamline has terminated.
   ind         = ind | all(~isfinite(T), 2);
   pos(~ind,1) = n(t(~ind));
   dir (ind)   = nan;

   % Change local coordinate when moving to new cell
   k = sub2ind(size(d), (1:size(dir,1))', dir);
   k = k(~isnan(k));
   pos(numel(dir) + k ) = 2-d(k);
end

% =========================================================================

function t = computeTime(xi, v)
% Compute time needed to reach xi=0 or xi=1 given velocities v=[v1,v2] at
% xi=0 and xi=1.  The formula is
%
%   t = xi/ui  or t = (1-xi)/ui,    if v1 = v2 = ui, and
%
%   t = 1/(v2-v1)*log(ue/ui),       otherwise
%
% where ui=v2*xi+v1*(1-xi) is the velocity at xi, and ue=v2 if ui>0 or
% ue=v1 if ui<0.

   tolerance = 100*eps;

   ui         = v(:,1) + xi.*diff(v, 1, 2);%(:,2)-v(:,1));
   ue         = v(:,    2);
   ue (ui<0)  = v(ui<0, 1);
   arg        = ue./ui;
   t          = inf(size(xi));

   % Linear velocity
   ind        = abs(diff(v, 1, 2)) > tolerance*abs(v(:,1));
   t(ind,:)   = 1./diff(v(ind,:), 1, 2).*log(arg(ind,:));

   % Constant velocity
   ds         = -xi;
   ds(ui > 0) = 1-xi(ui>0);
   t(~ind)    = ds(~ind)./ui(~ind);

   % nan happens for ui=ui=0
   t(arg<=0 | isnan(arg))   = inf;
end

% =========================================================================

function [x, i] = computePosition(xi, v, t)
% Compute position at time t given start point xi and velocities v=[v1,v2].
%
%   x = xi + v*t,    if v is constant or
%
%   x = xi + (ui*exp((v2-v1)*t) - ui)/(v2-v1), otherwise
%
   tolerance = 100*eps;

   du        = diff(v, 1, 2);
   ui        = v(:,1) + xi.*du;
   i         = 1 + ~(ui<0);
   x         = inf(size(xi));

   ind       = abs(du) > tolerance*abs(v(:,1));

   % linear velocity
   x(ind)    = xi(ind) + ( ui(ind).*exp(du(ind).*t(ind)) - ui(ind))./du(ind);

   % Constant velocity
   x(~ind,:) = xi(~ind,:) + v(~ind, 1).*t(~ind, :);
   x(~ind & t==inf) = xi(~ind & t==inf);
end
