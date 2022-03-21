function [G, msh] = gmshToMRST(filename)
% Create MRST grid object out of m file from Gmsh.
%
% SYNOPSIS:
%   G = gmshToMRST(filename)
%   [G, msh] = gmshToMRST(filename)
%
% DESCRIPTION:
%   Convert a Gmsh grid in the m-file filename to an MRST grid. The m file
%   from Gmsh can be generated in the Python interface by
%   gmsh.write("grid.m") or the GUI by saving it as a file with extension
%   .m.
%
% PARAMETERS:
%   filename - The .m file from Gmsh.
%
% RETURNS:
%   G   - MRST grid structure.
%   msh - Original Gmsh grid object.
%
% EXAMPLE:
%   Assuming there is a grid in Gmsh format in the file `grid.m`, this
%   grid can be converted to an MRST grid `G` by
%
%     G = gmshToMRST('grid.m');
%
%   The original Gmsh grid can be returned as an object named `msh` by
%
%     [G, msh] = gmshToMRST('grid.m');
%
% SEE ALSO:
%   `grid_structure`.

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

    run(filename);

    assert(exist('msh', 'var'), 'The file %s must construct a variable named msh', filename);

    if mrstVerbose
        disp(msh)
    end

    if msh.MIN(3) == msh.MAX(3)
        G = construct_2d(msh);
    else
        G = construct_3d(msh);
    end
    G.type = 'gmsh';

    G = computeGeometry(G, 'findNeighbors', true);
end


function G = construct_2d(msh)

    G.griddim = 2;

    % Set nodes
    G.nodes.num = msh.nbNod;
    G.nodes.coords = msh.POS(:, 1:G.griddim);

    % Cells are either triangles or quads
    cell_types = {'TRIANGLES', 'QUADS'};

    % Set up some constants (from Gmsh manual). All faces (edges) are
    % composed of 2 nodes
    edges = {[[0, 1]; [1, 2]; [2, 0]] + 1,...
             [[0, 1]; [1, 2]; [2, 3]; [3, 0]] + 1 };

    % Create a numbering of all the edges
    edgetable = [];
    for k = 1:numel(cell_types)
        if isfield(msh, cell_types{k})
            cell_nodes = msh.(cell_types{k});
            for e = 1:size(edges{k}, 1)
                edgetable = [edgetable; cell_nodes(:, edges{k}(e, 1:2))];
            end
        end
    end
    edgetable = unique(sort(edgetable, 2), 'rows');
    G.faces.num = size(edgetable, 1);

    % We need graph for G.cells.faces
    g = graph(edgetable(:, 1), edgetable(:, 2));

    % Construct the mrst data structures
    G.faces.nodes = reshape(edgetable, G.faces.num, 2).';
    G.faces.nodes = G.faces.nodes(:);
    G.faces.nodePos = (1:2:2*G.faces.num + 1)';

    % Initialize cells struct
    cells_num = zeros(numel(cell_types), 1);
    G.cells.faces = [];
    G.cells.facePos = 1;

    for k = 1:numel(cell_types)
        if isfield(msh, cell_types{k})
            cell_nodes = msh.(cell_types{k});
            num_local_faces = size(edges{k}, 1);
            edgenumbers = {};
            for e = 1:num_local_faces
                e1 = edges{k}(e, 1);
                e2 = edges{k}(e, 2);
                n1 = cell_nodes(:, e1);
                n2 = cell_nodes(:, e2);
                edgenumbers{e} = findedge(g, n1, n2);
            end
            edgenumbers = cell2mat(edgenumbers);
            cells_num(k) = size(edgenumbers, 1);

            % Construct the mrst data structures
            cells_faces = reshape(edgenumbers, cells_num(k), num_local_faces)';
            cells_faces = cells_faces(:);
            G.cells.faces = [G.cells.faces; cells_faces];

            cells_facePos = (num_local_faces:num_local_faces:num_local_faces*cells_num(k))';
            cells_facePos = cells_facePos + G.cells.facePos(end);
            G.cells.facePos = [G.cells.facePos; cells_facePos];
        end
    end
    G.cells.num = sum(cells_num);

    % Tags on lines
    G.faces.tag = zeros(G.faces.num, 1);
    if isfield(msh, 'LINES')
        edges = findedge(g, msh.LINES(:, 1), msh.LINES(:, 2));
        ii = find(edges > 0);
        G.faces.tag(edges(ii)) = msh.LINES(ii, 3);
    end

    % Tags on cells
    G.cells.tag = set_cell_tags(G, msh, cell_types, cells_num);

end


function G = construct_3d(msh)
    G.griddim = 3;

    % Set nodes
    G.nodes.num = msh.nbNod;
    G.nodes.coords = msh.POS(:, 1:G.griddim);

    % Cells are either tets, hexes, prisms or pyramids (the latter is
    % not implemented, because it has not yet been enountered)
    cell_types = {'TETS', 'HEXAS', 'PRISMS'};

    % Faces are either triangles or quads
    face_types = {'TRIANGLES', 'QUADS'};

    % Cell faces: tris and quads (exterior normal?)
    cells_faces_nodes = {
        { {[0, 1, 3], [0, 2, 1], [0, 3, 2], [1, 2, 3]}, {} },...
        { {}, {[0, 4, 7, 3], [1, 2, 6, 5], [0, 1, 5, 4], [2, 3, 7, 6], [0, 3, 2, 1], [4, 5, 6, 7]} },...
        { {[0, 2, 1], [3, 4, 5]}, {[0, 1, 4, 3], [0, 3, 5, 2], [1, 2, 5, 4]}}
                        };

    num_cell_faces = zeros(1, numel(cell_types));
    num_face_nodes = zeros(1, numel(face_types));
    for k = 1:numel(cell_types)
        for j = 1:numel(face_types)
            num_cell_faces(k) = num_cell_faces(k) + numel(cells_faces_nodes{k}{j});
            if ~isempty(cells_faces_nodes{k}{j})
                num_face_nodes(j) = numel(cells_faces_nodes{k}{j}{1});
            end
        end
    end
    assert(all(num_cell_faces == [4, 6, 5]));
    assert(all(num_face_nodes == [3, 4]));

    % Create numbering of the triangle and quad faces
    nodes_face = cell(2, 1);
    face_nodes = cell(2, 1);
    for j = 1:numel(face_types)
        % Create map from array of the face nodes (1 x 3 or 1 x 4) to
        % an integer number for the face.
        nodes_face{j} = containers.Map('KeyType', 'char', 'ValueType', 'uint64');
    end
    face_cnt = 0;

    for k = 1:numel(cell_types)
        if isfield(msh, cell_types{k})
            cell_nodes = msh.(cell_types{k});

            for j = 1:numel(face_types)
                % How many faces of this face type on this cell_type?
                num_local_faces = numel(cells_faces_nodes{k}{j});

                for i = 1:num_local_faces
                    local_nodes = cells_faces_nodes{k}{j}{i} + 1;
                    fn = cell_nodes(:, local_nodes);

                    if mrstVerbose
                        {cell_types{k}, face_types{j}, local_nodes}
                    end

                    for m = 1:size(fn, 1)
                        key = create_key(fn(m, :));

                        if ~isKey(nodes_face{j}, key)
                            face_cnt = face_cnt + 1;
                            dispif(mrstVerbose, ['new key ', key, ' of face nodes [', num2str(sort(fn(m, :))), '] face_cnt ', num2str(face_cnt), '\n']);
                            nodes_face{j}(key) = face_cnt;
                            face_nodes{j}(face_cnt, :) = cell_nodes(m, local_nodes);
                        elseif mrstVerbose
                            face_no = nodes_face{j}(key);
                            dispif(mrstVerbose, ['old key ', key, ' of face nodes [', num2str(face_nodes{j}(face_no, :)), '] face_no ', num2str(face_no), '\n']);
                        end

                    end

                end
            end
        end
    end

    G.faces.num = face_cnt;

    % Check that we have all found all faces
    nf = cell(numel(face_types), 1);
    for j = 1:numel(face_types)
        nf{j} = cell2mat(values(nodes_face{j}));
    end
    assert(all(sort(horzcat(nf{:})) == 1:G.faces.num));

    % Get faces to nodes
    G.faces.nodes = [];
    nodePos_diff = [];

    for j = 1:numel(face_types)
        fns = sort(cell2mat(values(nodes_face{j})));
        nodePos_diff = [nodePos_diff; repmat(num_face_nodes(j), numel(fns), 1)];
        for i = 1:numel(fns)
            fn = face_nodes{j}(fns(i), :);
            G.faces.nodes = [G.faces.nodes; fn.'];
        end
    end
    G.faces.nodePos = cumsum([1; nodePos_diff]);

    % Get cells to faces
    cells_num = zeros(numel(cell_types), 1);
    G.cells.faces = [];
    facePos_diff = [];

    for k = 1:numel(cell_types)
        if isfield(msh, cell_types{k})
            cell_nodes = msh.(cell_types{k});
            cells_num(k) = size(cell_nodes, 1);
            facePos_diff = [facePos_diff; repmat(num_cell_faces(k), size(cell_nodes, 1), 1)];

            for ci = 1:cells_num(k)
                for j = 1:numel(face_types)
                    num_local_faces = numel(cells_faces_nodes{k}{j});

                    for i = 1:num_local_faces
                        nodes = cells_faces_nodes{k}{j}{i} + 1;
                        fn = cell_nodes(ci, nodes);
                        key = create_key(fn);
                        face_no = nodes_face{j}(key);
                        dispif(mrstVerbose, [key, ' ', num2str(fn), ' ', num2str(face_no), '\n'])
                        G.cells.faces = [G.cells.faces; face_no];
                    end
                end
            end
        end
    end

    G.cells.num = sum(cells_num);
    G.cells.facePos = cumsum([1; facePos_diff]);

    % Tags on triangles or quads
    G.faces.tag = zeros(G.faces.num, 1);
    for j = 1:numel(face_types)
        if isfield(msh, face_types{j})
            fns = msh.(face_types{j});

            % Make sure we have tags
            assert(size(fns, 2) == num_face_nodes(j) + 1);

            for i = 1:size(fns, 1)
                key = create_key(fns(i, 1:end - 1));
                face_no = nodes_face{j}(key);
                G.faces.tag(face_no) = fns(i, end);
            end
        end
    end

    % Tags on cells
    G.cells.tag = set_cell_tags(G, msh, cell_types, cells_num);
    
end


function key = create_key(fn)
    fn = sort(fn);
    key = md5sum(fn);
end


function cell_tags = set_cell_tags(G, msh, cell_types, cells_num)

    cell_tags = zeros(G.cells.num, 1);
    cidx = [0; cumsum(cells_num)];
    for k = 1:numel(cell_types)
        if isfield(msh, cell_types{k})
            cell_tags((cidx(k) + 1):cidx(k + 1)) = msh.(cell_types{k})(:, end);
        end
    end

end
