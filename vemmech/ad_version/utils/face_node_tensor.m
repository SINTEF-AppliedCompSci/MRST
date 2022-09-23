function tensor = face_node_tensor(G, fname, nname, varargin)
% Construct the face-node indicator tensor for grid G 
% 
% SYNOPSIS:
%   function tensor = face_node_tensor(G)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G     - grid structure
%   fname - name to use for face index (e.g. 'f')
%   nname - name to use for node index (e.g. 'n')
% 
% OPTIONAL PARAMETERS:
%   'values'        - vector with values to use instead of 1 to indicate each 
%                     face-node coincidence.  Can be used to directly
%                     construct an arbitrary tensor in face-node space.  
%                     NB: It is the user's responsibility to ensure the
%                     vector has the correct number of elements.
%   'boundary_only' - eliminate all non-boundary faces.
%
% RETURNS:
%   tensor - indicator tensor with value '1' for all combinations of faces
%            and nodes that are neighbors, and '0' for all other combinations
%
% EXAMPLE:
%  fnt = face_node_tensor(G, 'face', 'node')

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

opt = struct('values', [], 'boundary_only', false);
opt = merge_options(opt, varargin{:});

ixs = [rldecode((1:G.faces.num)', diff(G.faces.nodePos)), G.faces.nodes];

if opt.boundary_only
   
   % faces to keep
   int_faces = prod(G.faces.neighbors,2) == 0;
   
   keep = rldecode(int_faces, diff(G.faces.nodePos));
   ixs = ixs(keep,:);
   
   if ~isempty(opt.values)
      opt.values = opt.values(keep);
   end
   
end

tensor = SparseMultiArray(opt.values, ixs, {fname, nname});
