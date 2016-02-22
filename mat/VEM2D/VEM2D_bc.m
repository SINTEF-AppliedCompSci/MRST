function [bcDof, bBC] = VEM2D_bc(G, bc, k)
%--------------------------------------------------------------------------
%   Sets boundary conditions for the Poisson problem.
%
%   G:  MRST grid.
%   bc: Struct of boundary conditions.
%       filed1: 'bcFunc', boundary condition function handles.
%       filed2: 'bcFaces', vectors of boundary faces.
%       fild3:  'bcType', 'neu' for Neumann, 'dir' for Dirichlet.
%       Example:
%
%       bc = struct('bcFunc', {{gN, gD}},'bcFaces', ...
%           {{boundaryEdges(bNeu), boundaryEdges(~bNeu)}}, ...
%           'bcType', {{'neu', 'dir'}});
%--------------------------------------------------------------------------

nN = G.nodes.num;           %   Number of nodes.
nE = G.faces.num;           %   Number of edges.
nK = G.cells.num;           %   Number of cells.

if k == 1
    NK = nN;
elseif k == 2
    NK = nN + nE + nK;
end

bcFunc = bc.bcFunc;         %   Boundary condition functions.
bcFaces = bc.bcFaces;       %   Vectors of boundary edges.
bcType = bc.bcType;         %   Boundary condition types
Nbc = numel(bcFunc);        %   Number of boundary conditions.
bcDof = zeros(NK, 1);     %   Dof map for boundary conditions.
bBC = zeros(NK, 1);       %   b vector for boundary conditions.
for b = 1:Nbc
    g = bcFunc{b};          %   BC funciton.
    edges = bcFaces{b};     %   BC edges.
    if size(edges,1) == 1
        edges = edges';
    end
    type = bcType{b};       %   BC type.
    n = numel(edges);       %   Number of edges where BC applies.
        
    nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges+1)-1);
    nodes = G.faces.nodes(nodeNum);
    if size(nodes,1) == 1
        nodes = nodes';
    end
    if k == 1
        X = G.nodes.coords(nodes,:);
        dofVec = nodes';
    elseif k == 2                     %   Dof coordinates
        X = [G.nodes.coords(nodes,:); G.faces.centroids(edges,:)];
        dofVec = [nodes', edges' + nN];
    end
                            %   Fix dimensions.
    if strcmp(type, 'dir')
                            %   Apply Dirichelt BC's.
        bcDof(dofVec) = 1;
        bBC(dofVec) = g(X);
    elseif strcmp(type, 'neu')
                            %   Apply Neuman BC's.
        bcDof(dofVec) = 2;
        edgeLengths = G.faces.areas(edges);
                            %   Contribution from each edge to
                            %   \int_\partial \Omega g_N \phi_i ds
        for e = 1:n
            bBC([nodes(2*e-1), nodes(2*e), nN + edges(e)]) = ...
                bBC([nodes(2*e-1), nodes(2*e), nN + edges(e)]) + ...
                edgeLengths(e).*[1/6*g(X(2*e-1,:)); 1/6*g(X(2*e,:)); 2/3*g(X(2*n + e,:))];
        end
    end
end

end