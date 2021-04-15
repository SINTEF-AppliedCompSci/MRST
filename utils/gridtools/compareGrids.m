function compareGrids(G1, G2)
%Determine if two grid structures are the same.
%
% SYNOPSIS:
%   compareGrids(G1, G2)
%
% DESCRIPTION:
%   This function considers two grids to be the same if all fundamental
%   grid_structure fields have the same size and values, and in the same
%   order within each field, in both grids.
%
%   Otherwise, if the grids contain equally sized and valued fundamental
%   fields, but the values differ in internal ordering, then the grids are
%   topologically equivalent but differently represented internally.
%
%   Otherwise, the grids are not equivalent in any sense.
%
% PARAMETERS:
%   G1,G2 - Two grid structures as defined in 'grid_structure'.
%
% RETURNS:
%   Nothing - Diagnostic messages are printed to standard output stream.
%
% NOTE:
%   Function `compareGrids` is mainly a grid processing debugging tool.
%
% SEE ALSO:
%   `grid_structure`.

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


   hash1 = md5sum(G1);
   hash2 = md5sum(G2);
   if all(hash1 == hash2),
      fprintf('Grids are (binary) identical\n');
      return;
   end

   tol = 10 * eps(max(abs([G1.nodes.coords(:); G2.nodes.coords(:)])));

   a = compareFields(G1, G2, tol);
   if ~all(a),

      G1 = sortGrid(G1);
      G2 = sortGrid(G2);

      if any([G1.nodes.num ~= G2.nodes.num, ...
              G1.faces.num ~= G2.faces.num, ...
              G1.cells.num ~= G2.cells.num]),

          fprintf('Grids are not equivalent.\n');
          if G1.nodes.num ~= G2.nodes.num,
              fprintf(' * The number of nodes differ.\n');
          end
          if G1.faces.num ~= G2.faces.num,
              fprintf(' * The number of faces differ.\n');
          end
          if G1.cells.num ~= G2.cells.num,
              fprintf(' * The number of cells differ.\n');
          end

      else

         a = compareFields(G1, G2, tol);
         if all(a),
            fprintf(['Grids have different representation, ', ...
                     'but are topologically equivalent.\n']);
         else
            fprintf('Grids are not equivalent because\n');
            if ~a(1), fprintf(' * nodes.coords differ\n');    end
            if ~a(2), fprintf(' * faces.nodePos differ\n');   end
            if ~a(3), fprintf(' * faces.neighbors differ\n'); end
            if ~a(4), fprintf(' * faces.nodes differ\n');     end
            if ~a(5), fprintf(' * cells.facePos differ\n');   end
            if ~a(6), fprintf(' * cells.indexMap differ\n');  end
            if ~a(7), fprintf(' * cells.faces differ\n');     end
         end
      end

      % In a valid grid, any node in a cell should occur more that 2 times
      % to ensure there are no cell with holes.
   else
      fprintf('Grids are identical to tolerance %g in node coordinates.\n',...
         tol);
   end
end

%--------------------------------------------------------------------------

function a = compareFields(G1, G2, tol)
   a = [norm(G1.nodes.coords(:) - G2.nodes.coords(:))<tol,...
        G1.faces.num == G2.faces.num && all(G1.faces.nodePos      == G2.faces.nodePos),   ...
        G1.faces.num == G2.faces.num && all(G1.faces.neighbors(:) == G2.faces.neighbors(:)),   ...
        all(size(G1.faces.nodes)==size(G2.faces.nodes)) && all(G1.faces.nodes == G2.faces.nodes),             ...
        all(G1.cells.facePos      == G2.cells.facePos),...
        all(G1.cells.indexMap     == G2.cells.indexMap),...
        all(G1.cells.faces(:,1)   == G2.cells.faces(:,1))];
end
