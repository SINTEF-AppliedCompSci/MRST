function [X, Y, Z, lineIx] = buildCornerPtNodes(grdecl, varargin)
%Construct physical nodal coordinates for CP grid.
%
% SYNOPSIS:
%   [X, Y, Z, lineIx] = buildCornerPtNodes(grdecl)
%   [X, Y, Z, lineIx] = buildCornerPtNodes(grdecl, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl  - Eclipse file output structure as defined by `readGRDECL`.
%             Must contain at least the fields `cartDims`, `COORD` and
%             `ZCORN`.
%
% KEYWORD ARGUMENTS:
%
%  'Verbose'              - Whether or not to emit informational messages.
%                           Default value: Verbose = mrstVerbose.
%
%  'CoincidenceTolerance' - Absolute tolerance used to detect collapsed
%                           pillars where the top pillar point coincides
%                           with the bottom pillar point.  Such pillars are
%                           treated as is they were vertical.
%                           Default value: CoincidenceTolerance = `100*eps`.
%
% RETURNS:
%   X, Y, Z - Size `2*grdecl.cartDims` arrays of `x`, `y` and `z` physical
%             nodal coordinates, respectively.
%
%   lineIx  - Index of pillar line. There are (nx+1)x(ny+1) pillar lines
%             along which node coordinates are defined by their `ZCORN`
%             value.
% EXAMPLE:
%   gridfile  = [DATADIR, filesep, 'case.grdecl'];
%   grdecl    = readGRDECL(gridfile);
%   [X, Y, Z] = buildCornerPtNodes(grdecl);
%
% SEE ALSO:
%   `readGRDECL`, `processGRDECL`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('Verbose', mrstVerbose, 'CoincidenceTolerance', 100*eps);
   opt = merge_options(opt, varargin{:});

   % Recover logical (cell) dimension of grid
   nx = grdecl.cartDims(1);
   ny = grdecl.cartDims(2);
   nz = grdecl.cartDims(3);

   % Enumerate individual pillars in grid
   % Pillars 'p1', 'p2', 'p3' and 'p4' bound individual cells in grid.
   pillarIx = reshape(1 : (nx+1)*(ny+1), [nx+1, ny+1]);

   p1 = pillarIx(1:nx,   1:ny);
   p2 = pillarIx(2:nx+1, 1:ny);
   p3 = pillarIx(1:nx,   2:ny+1);
   p4 = pillarIx(2:nx+1, 2:ny+1);
   clear pillarIx

   shape = 2 .* [nx, ny, nz];
   lineIx = zeros(shape);
   lineIx(1:2:2*nx, 1:2:2*ny, :) = repmat(p1, [1, 1, 2*nz]);
   lineIx(2:2:2*nx, 1:2:2*ny, :) = repmat(p2, [1, 1, 2*nz]);
   lineIx(1:2:2*nx, 2:2:2*ny, :) = repmat(p3, [1, 1, 2*nz]);
   lineIx(2:2:2*nx, 2:2:2*ny, :) = repmat(p4, [1, 1, 2*nz]);
   clear p1 p2 p3 p4

   lines = reshape(grdecl.COORD, 6, []) .';
   lines = lines(lineIx(:), :);
   %clear lineIx

   % Recover physical nodal coordinates along pillars
   %
   % We assume that all pillars are straight lines so linear interpolation
   % is sufficient.
   %
   % As a special case we need to handle "collapsed" pillars of the form
   %
   %    x0 y0 z0   x0 y0 z0
   %
   % where the top pillar point coincides with the bottom pillar point.
   % Following ECLIPSE, we assume that such pillars are really vertical...
   %
   ix    = abs(lines(:,6) - lines(:,3)) < abs(opt.CoincidenceTolerance);
   t     = (grdecl.ZCORN(:) - lines(:,3)) ./ (lines(:,6) - lines(:,3));
   t(ix) = 0;

   nCollapsed = sum(ix);
   dispif(opt.Verbose && nCollapsed > 0, ...
         ['Detected %d collapsed (of %d total) pillars.\n', ...
          'Treated as vertical...\n'], ...
          nCollapsed, prod(grdecl.cartDims(1:2) + 1));

   xCoords = lines(:,1) + t.*(lines(:,4) - lines(:,1));
   yCoords = lines(:,2) + t.*(lines(:,5) - lines(:,2));
   clear t ix lines

   % Reshape to form more amenable to further processing
   X = reshape(xCoords,      shape); clear xCoords
   Y = reshape(yCoords,      shape); clear yCoords
   Z = reshape(grdecl.ZCORN, shape);
end
