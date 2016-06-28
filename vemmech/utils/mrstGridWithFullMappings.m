function G = mrstGridWithFullMappings(G)
% Add all mappings to mrst grid. This is needed for VEM of finite element
% type methods. 
%
% SYNOPSIS:
%   G = mrstGridWithFullMappings(G)
%   
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% OPTIONAL PARAMETERS:
%
%   'hingenodes'  - A struct with fields 'faces' and 'nodes'.  A hinge node
%                   is an extra center node for a face, that is used to
%                   triangulate the face geometry.  For each face number F
%                   in 'faces', there is a row in 'nodes' which holds the
%                   node coordinate for the hinge node belonging to face F.
%
%
% COMMENTS:
%   At the moment, mappings that are not needed may be constructed. The mapping
%   from local nodes to global is not calculated and needs to be made localy
%   if needed,  for example by using a sparse matrix.
%   The methods assume that edges and faces are sorted.
%
% SEE ALSO:
%   grid_structure. calulateGeometryCalc,  VEM_linElast
%

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 

    if(G.griddim == 3)
        ind1 = mcolon(G.faces.nodePos(1:end-1), G.faces.nodePos(2:end)-1, 1);
        ind2 = mcolon(G.faces.nodePos(1:end-1)+1, G.faces.nodePos(2:end), 1);
        ind2(G.faces.nodePos(2:end)-1) = G.faces.nodePos(1:end-1);
        eh = [G.faces.nodes(ind1), G.faces.nodes(ind2)];
        ehs = sort(eh, 2);
        [e, ~, j] = unique(ehs, 'rows', 'sorted');
        assert(all(all(e(j, :) == ehs)))
        
        G.faces.edges = j;
        G.faces.edgePos = [1;cumsum(diff(G.faces.nodePos))+1];
        G.faces.edgeSign = 2*(e(j, 1) == eh(:, 1))-1;
        
        % make ting si pos format
        e = e';
        ne = size(e', 1);
        G.edges.nodes = e(:);
        nPos = cumsum(repmat(2, ne, 1));
        G.edges.nodePos = [1;nPos+1];
        G.edges.num = ne;
    else
        faces = G.cells.faces(:, 1);
        inodes = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1);
        assert(2*numel(faces) == numel(inodes))
        cell = rldecode([1:G.cells.num]', diff(G.cells.facePos));%#ok
        sign = 2*(cell == G.faces.neighbors(faces, 1))-1;
        nf = [1:2:numel(inodes)]'-(sign-1)/2;
        nodes = G.faces.nodes(inodes(nf));
        G.cells.nodes = nodes;
        G.cells.nodePos = G.cells.facePos;
    end
    w = createGridMappings(G);


    %% Setup node to face mapping
    n2f = unique([w(:, 2), w(:, 4)], 'rows', 'sorted');
    [pos, val] = map2Pos(n2f);
    G.nodes.facePos = pos;
    G.nodes.faces = val;


    %% Setup cell to node mapping
    if(G.griddim == 3)
        c2n = unique([w(:, 1), w(:, 2)], 'rows', 'sorted');
        [pos, val] = map2Pos(c2n);
        G.cells.nodePos = pos;
        G.cells.nodes = val;
    end

    %% Setup node to cell mapping
    n2c = unique([w(:, 2), w(:, 1)], 'rows', 'sorted');
    [pos, val] = map2Pos(n2c);
    G.nodes.cellPos = pos;
    G.nodes.cells = val;
    G.type = [G.type, { mfilename }];

end

function [pos, val] = map2Pos(nn)
    ind = find(diff(nn(:, 1)) == 1);
    pos = [1;ind+1;size(nn, 1)+1];
    val = nn(:, 2);
end

