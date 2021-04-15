function G = mcomputeGeometry(G)
%Compute geometric primitves using compiled C code.
%
% SYNOPSIS:
%   G = mcomputeGeometry(G)
%
% PARAMETERS:
%   G - A grid_structure.
%
% RETURNS:
%   G - An updated grid_structure containing areas, volumes, normals and
%       centroids.
%
% SEE ALSO:
%   `computeGeometry`, `grid_structure`.

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


   if numel(G) > 0
      if ~isfield(G(1), 'type')
         warning(msgid('GridType:Unknown'),                         ...
                ['Input grid has no known type. ',                  ...
                 'I''ll assume it arose from the primordial soup...']);

         [ G.type ] = deal( {'Primordial Soup'} );
      end

      for k = 1 : numel(G)
         if size(G(k).nodes.coords, 2) ~= 3
            error(['Function ''%s'' is only supported in three ', ...
                   'space dimensions.'], mfilename);
         end

         [fa, fc, fn, cc, cv] = mex_compute_geometry(G(k));

         if ~ all(fa > 0)
            warning(msgid('FaceAre:NonPositive'), ...
                    'Face area negative or zero');
         end
         if ~ all(cv > 0)
            warning(msgid('CellVolume:NonPositive'), ...
                    'Cell volume negative or zero');
         end

         G(k).faces.areas     = fa;
         G(k).faces.centroids = fc .';
         G(k).faces.normals   = fn .';
         G(k).cells.centroids = cc .';
         G(k).cells.volumes   = cv;

         G(k).type = [ G(k).type, { mfilename } ];
      end
   end
end
