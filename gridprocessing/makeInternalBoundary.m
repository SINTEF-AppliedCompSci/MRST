function [G, N] = makeInternalBoundary(G, faces, varargin)
%Make internal boundary in grid along FACES
%
% SYNOPSIS:
%   G = makeInternalBoundary(G, f)
%   G = makeInternalBoundary(G, f, 'pn', pv)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   f       - (unique) faces along which a boundary will be inserted.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   tag    -  Tag inserted in G.faces.tag for faces on the internal
%             boundary.
%
% RETURNS:
%   G       - Modified grid structure.
%
%   N       - An n x 2 array of face numbers.  Each pair in the array
%             correspond to a face in FACES.  N is sufficient to create
%             flow over the internal boundary or to remove the boundary.
%
% COMMENTS:
%
% SEE ALSO:
%  removeInteralBoundary

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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



opt = struct('tag', -1);
opt = merge_options(opt, varargin{:});

ind = all(G.faces.neighbors(faces,:)~=0, 2);
n   = numel(faces(ind));

if n == 0,
    error('No internal faces given')
end

% Repeated faces generate exotic grid errors.
uf = unique(faces);
if numel(faces) ~= numel(uf),
   warning('List of faces contains repeated face indices.');
   faces = uf;
   clear uf
end

% Copy face info
[fnodes, nnodes, neigh, tags] = copyFaces(G, faces(ind));

% Remove faces
G      = removeFaces(G, faces(ind));

% Double number of faces, introducing 0-neighbors
neigh  = reshape([neigh(:,1), zeros(n,2), neigh(:,2)]', 2, [])';
if ~isempty(tags),
   tags   = reshape([tags(:,1), zeros(n,2), tags(:,2)]', 2, [])';
end

pos    = cumsum([1;double(nnodes)]);
p1     = rldecode(pos(1:end-1), 2*ones(n,1));
p2     = rldecode(pos(2:end)-1, 2*ones(n,1));
ix     = mcolon(p1,p2);
nnodes = rldecode(nnodes, 2*ones(n,1));
fnodes = fnodes(ix);

% Add replacement faces
N      = G.faces.num+reshape(1:numel(neigh)/2, 2, []) .';
if any(tags)
   [G,new_faces]      = addFaces(G, fnodes, nnodes, neigh, tags);
   G.faces.tag(new_faces) = opt.tag;
else
   [G,new_faces]      = addFaces(G, fnodes, nnodes, neigh);
end

G.type = [G.type, { mfilename }];
