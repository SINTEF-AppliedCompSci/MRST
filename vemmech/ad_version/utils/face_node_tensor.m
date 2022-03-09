function tensor = face_node_tensor(G, fname, nname, varargin)
%
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
%
% SEE ALSO:
%
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
