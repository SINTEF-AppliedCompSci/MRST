function [G, N] = makeInternalBoundary(G, faces, varargin)
%Make internal boundary in grid along specified faces.
%
% SYNOPSIS:
%    H     = makeInternalBoundary(G, faces)
%    H     = makeInternalBoundary(G, faces, 'pn1', pv1, ...)
%   [H, N] = makeInternalBoundary(...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   faces   - Faces along which a boundary will be inserted.  Vector of
%             face indices.  Repeated indices are supported but generally
%             discouraged.  One pair of new grid faces will be created in
%             the output grid for each unique face in 'faces' provided the
%             identified face is not located on the boundary of 'G'.
%
% KEYWORD ARGUMENTS:
%
%  'tag'   - Tag inserted into 'G.faces.tag' (if present) for faces on the
%            internal boundary.  Integer scalar. 
%            Default value: tag = -1.
%
%  'check' - Whether or not to check for repeated face indices in input
%            `faces`.  Emit diagnostic when set and in presence of repeated
%            indices.  Logical scalar.
%            Default value: check = LOGICAL(mrstVerbose).
%
% RETURNS:
%   H - Modified grid structure.
%
%   N - An n-by-2 array of face indices. Each pair in the array corresponds
%       to a face in `faces`.  In particular, `faces(i)` in the input grid
%       `G` is replaced by the pair of faces `[N(i,1), N(i,2)]` in the
%       result grid `H`.
%
%       If any of the faces in `faces` are on the boundary, then the
%       corresponding rows in `N` is `NaN` in both columns.
%
%       The quantity `N` is sufficient to create flow over the internal
%       boundary or to remove the internal boundary at some later point.
%
% SEE ALSO:
%   `removeInternalBoundary`, `mrstVerbose`

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

opt = struct('tag', -1, 'check', logical(mrstVerbose));
opt = merge_options(opt, varargin{:});

[faces, iu] = unique_faces(faces, inputname(2), opt);

int = ~ any(G.faces.neighbors(faces,:) == 0, 2);
n   = sum(int);
if n == 0,
   % All of the supplied 'faces' are on the boundary.  Don't change input
   % grid and return an all-NaN 'N' result.
   warning(msgid('IFaces:NotSupplied'), 'No internal faces supplied.');
   N = NaN([numel(iu), 2]);
   return
end

% Extract face info for user-specified, internal faces.
faces = faces(int);
[fnodes, nnodes, neigh, tags] = copyFaces(G, faces);

% Double number of faces, introducing zero (outside) neighbours/tags.
expand = @(a) reshape([a(:,1), zeros([n, 2]), a(:,2)] .', 2, []) .';
neigh  = expand(neigh);
if ~isempty(tags),
   tags = expand(tags);
end

% Remove user-specified faces.
G = removeFaces(G, faces);



pos    = cumsum([1;double(nnodes)]);
p1     = rldecode(pos(1:end-1), 2*ones(n,1));
p2     = rldecode(pos(2:end)-1, 2*ones(n,1));
ix     = mcolon(p1,p2);
nnodes = rldecode(nnodes, 2*ones(n,1));
fnodes = fnodes(ix);
% Add replacement faces.
N         = NaN([numel(int), 2]);
N(int, :) = G.faces.num + reshape(1 : (numel(neigh) / 2), 2, []) .';
N         = N(iu, :);

if any(tags),
   [G, new_faces]         = addFaces(G, fnodes, nnodes, neigh, tags);
   G.faces.tag(new_faces) = opt.tag;
else
   G = addFaces(G, fnodes, nnodes, neigh);
end

G.type = [G.type, { mfilename }];
end

%--------------------------------------------------------------------------

function [faces, iu] = unique_faces(faces, iname, opt)
   select_1st_nonempty_str = @(varargin) ...
      varargin { ...
      find(cellfun(@(x) ~isempty(x) && ischar(x), varargin), 1) };

  % Repeated faces generate exotic grid errors.
  [faces, iu, iu] = unique(faces);                              %#ok<ASGLU>

  dispif(opt.check && (numel(faces) < numel(iu)), ...
         ' * Input ''%s'' contains repeated indices. *\n', ...
         select_1st_nonempty_str(iname, 'faces'));
end
