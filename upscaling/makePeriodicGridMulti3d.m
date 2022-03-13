function [Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dprl, varargin)
% Construct a periodic grid.
% G grid, bcl left boundary conditions, bcr right boundary conditions, dprl
% delta pressure over periodic edges

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

% New cells.faces will have an additional column
faces = G.cells.faces(:, [1, 1, 2:end]); % Duplicate column 1 to allocate.

% Insert tags for new periodic boundary at position 2
faces(:, 2) = (1 : size(faces, 1)) .';

G.cells.faces = faces;

N_bcl = numel(bcl);
% Assert that the numbers of boundary conditions are all equal.
assert(isequal(N_bcl, numel(bcr), numel(dprl)));

rm_faces = cell ([N_bcl, 1]);
cells    = cell ([N_bcl, 1]);
ns       = zeros([N_bcl, 1]);

G_fc = G.faces.centroids;

for i = 1 : N_bcl
    % We want to iterate over principal directions which are not the
    % current direction.
    d = 1:G.griddim;
    d = d(d ~= i);
    assert(i <= G.griddim)

    % Match up indexes by sorting the centroids of the active principal
    % directions
    [~, ia] = sortrows(G_fc(bcl{i}.face, d));
    [~, ib] = sortrows(G_fc(bcr{i}.face, d));

    faces1 = bcl{i}.face(ia);
    faces2 = bcr{i}.face(ib);

    if G.griddim > 1
        % Verify that we have identified the same faces. This is done by
        % checking centroids and face areas.
        assert(all(sum((G_fc(faces1,d) - G_fc(faces2,d)).^2, 2) < 1e-10))
        assert(all(abs(G.faces.areas(faces1)-G.faces.areas(faces2))<1e-10))
    end

    rm_faces{i} = [faces1, faces2];
    ns(i) = numel(faces1);
    cells{i} = sum(G.faces.neighbors(faces1,:),2);
end

dp = vertcat(dprl{:});

% Remove the boundary between the boundaries. rm_faces when expanded will
% be a nx2 vector where faces in the first column will be connected to
% faces in the second column.
[Gp, f] = removeInternalBoundary(G, vertcat(rm_faces{:}));
cells = vertcat(cells{:});

dp = rldecode(dp, ns);

tags = rldecode((1:G.griddim)', ns);
bcp  = struct('face',    f,    ...
              'value',   dp,   ...
              'tags',    tags, ...
              'sign',    2*(Gp.faces.neighbors(f,1) == cells) - 1, ...
              'upcells', cells);

Gp.cells.faces = [Gp.cells.faces(:,1), ...
                  G.cells.faces(Gp.cells.faces(:,2), 3:end), ...
                  Gp.cells.faces(:,2)];

Gp.cells.centroids = G.cells.centroids;
Gp.faces.normals   = [];
Gp.nodes.coords    = [];

%assert(~(any(any(Gp.faces.neighbors==0,2))));
% remap areas of faces assume periodic boundary is made on faces with
% the same area.
Gp.faces.areas(Gp.cells.faces(:,1)) = ...
   G.faces.areas(G.cells.faces(Gp.cells.faces(:,end), 1));

% Helper function 'removeInternalBoundary' does not touch Gp.faces.areas.
% Ensure that we have correctly sized arrays in the output grid.
Gp.faces.areas = Gp.faces.areas(1 : Gp.faces.num);
end
