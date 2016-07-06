function [div] = VEM3D_div(G)
%
%
% SYNOPSIS:
%   function [div] = VEM3D_div(G)
%
% DESCRIPTION:
% Computes a discrete div operator in 3D. This discrete div operator is a
% mapping from node-valued displacement vector to cell-valued 2D
% vector. Node-valued displacement vectors correspond to the degrees of freedom
% that determine for each cell a displacement function over the cell via the
% virtual basis function. The discrete div operator computes the L^2
% projection, cell-wise,  of this displacement function. For more detail, see
% ecmor paper.
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   div - matrix corresponding to the discrete div operator.
%
% EXAMPLE:
%
% SEE ALSO:
%


    assert(G.griddim == 3);
    cells     = 1:G.cells.num;
    ifaces    = mcolon(G.cells.facePos(cells), (G.cells.facePos(cells + 1) - 1));
    faces     = G.cells.faces(ifaces, 1);
    inodes    = mcolon(G.faces.nodePos(faces), (G.faces.nodePos(faces + 1) - 1));
    nodes     = G.faces.nodes(inodes);
    cellfaces = rldecode(cells', diff(G.cells.facePos));
    nvec      = bsxfun(@times, G.faces.normals(faces, :), ...
                       (2 * (G.faces.neighbors(faces, 1) == cellfaces) - 1));
    nfl       = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces);
    nvec  = bsxfun(@rdivide, nvec, G.faces.areas(faces));
    [qc, qf] = calculateQC(G);
    nvec = bsxfun(@rdivide, nvec, G.faces.areas(faces));
    matel = reshape(bsxfun(@times, rldecode(nvec, nfl), qf(inodes))', [], 1);
    jind  = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim);
    iind  = rldecode(cellfaces, nfl * 3);
    div   = sparse( iind', jind', matel(:), G.cells.num, G.nodes.num .* G.griddim);

end
