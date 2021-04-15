function h = plotFaceOutline(G, f, varargin)
%PLOTFACEOUTLINE plots the outline of a selection of grid faces.
%
% SYNOPSIS:
%       plotFaceOutline(G, faces)
%       plotFaceOutline(G, faces, 'pn1', pv1, ...)
%   h = plotFaceOutline(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   faces   - Vector of face indices.  The graphical output of
%             'plotFaceOutline' will be restricted to the subset of grid
%             faces from 'G' represented by 'faces'.
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PLOT3
%             meaning all properties supported by PLOT3 are valid.
%
% RETURNS:
%   h - Handle to resulting patch objects.  The patch objects are added to
%       the current AXES object.
%
% NOTES:
%   Function 'plotFaceOutline' is implemented directly in terms of the
%   basic plotting function PATCH to enable control via the same set of
%   options as in previous implementations of the plotGrid family.
%
% EXAMPLE:
%   % Plot grid with boundary faces on left side in red colour:
%   G     = cartGrid([5, 5, 2]);
%   faces = boundaryFaceIndices(G, 'LEFT', 1:5, 1:2, []);
%   plotGrid (G, 'faceColor', 'none'); view(3)
%   plotFaceOutline(G, faces, 'r');
%
% SEE ALSO:
%   `plotCellData`, `plotGrid`, `plotFaces`, `newplot`, `patch`.


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



   % For each collection of faces f(pos(i):pos(i+1)-1), plot edges on
   % boundary of collection
   if isfield(G, 'parent')
      if G.griddim == 3
         pos     = cumsum(double([1; G.faces.connPos(f+1)-G.faces.connPos(f)]));
         ix      = mcolon(G.faces.connPos(f), G.faces.connPos(f+1)-1);
         f       = G.faces.fconn(ix);
         cno     = rldecode(1:numel(pos)-1, diff(pos), 2)';
      else
         ix    = false(G.cells.num,1);
         ix(f) = true;
         f     = find(ix(G.partition));
         cno   = G.partition(f);
      end
      G       = G.parent;
   end

   if G.griddim == 3
      cno    = rldecode(cno, ...
                        G.faces.nodePos(f+1)-G.faces.nodePos(f));
      ix     = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
      e(:,1) = G.faces.nodes(ix);
      e(:,2) = e([2:end, 1],1);

      pos    = cumsum(double([1; G.faces.nodePos(f+1) - G.faces.nodePos(f)]));
      e(pos(2:end)-1, 2) = e(pos(1:end-1), 1);

   elseif G.griddim ==2
      cno = rldecode(cno, G.cells.facePos(f+1)-G.cells.facePos(f));
      ix  = mcolon(G.cells.facePos(f), G.cells.facePos(f+1)-1);
      f   = G.cells.faces(ix);
      ix  = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
      e   = reshape(G.faces.nodes(ix), 2, [])';

   else
      error('Field ''griddim'' should be 2 or 3.');
   end

   % The standard phrase to build a new topological map...
   E      = [sort(e,2), cno];
   E      = sortrows(E, [3,1,2]);
   [~, n] = rlencode(E);
   N      = rldecode(n, n); 

   h      = plotLineSegments(G, E(N==1, 1:2), varargin{:});
end

function h = plotLineSegments(G, e, varargin)
% Plot all line segments given by node pairs in each row in e.
   e = unique([e;fliplr(e)], 'rows');
   h = patch('vertices', G.nodes.coords, 'faces', e, varargin{:});
end

