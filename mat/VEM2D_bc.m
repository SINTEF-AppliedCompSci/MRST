function [bcDof, bBC] = VEM2D_bc(G, bc)
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

Nn = G.nodes.num;           %   Number of nodes.
Ne = G.faces.num;           %   Number of edges.
Nc = G.cells.num;           %   Number of cells.

Ndof = Nn + Ne + Nc;        %   Number of dofs.

bcFunc = bc.bcFunc;         %   Boundary condition functions.
bcFaces = bc.bcFaces;       %   Vectors of boundary edges.
bcType = bc.bcType;         %   Boundary condition types
Nbc = numel(bcFunc);        %   Number of boundary conditions.
bcDof = zeros(Ndof, 1);     %   Dof map for boundary conditions.
bBC = zeros(Ndof, 1);       %   b vector for boundary conditions.
for b = 1:Nbc
    g = bcFunc{b};          %   BC funciton.
    edges = bcFaces{b};     %   BC edges.
    type = bcType{b};       %   BC type.
    n = numel(edges);       %   Number of edges where BC applies.
    nodes = zeros(2*n,1);   %   Nodes where BC applies.
    for e = 1:n
        nodeNum = G.faces.nodePos(edges(e)):G.faces.nodePos(edges(e) + 1)-1;
        nodes(2*e-1:2*e) = G.faces.nodes(nodeNum);
    end
                            %   Dof coordinates
    X = [G.nodes.coords(nodes,:); G.faces.centroids(edges,:)];
                            %   Fix dimensions.
    edges = fixDim(edges);
    nodes = fixDim(nodes);
    if strcmp(type, 'dir')
                            %   Apply Dirichelt BC's.
        bcDof([nodes, Nn + edges]) = 1;
        bBC([nodes, Nn + edges]) = g(X);
    elseif strcmp(type, 'neu')
                            %   Apply Neuman BC's.
        bcDof([nodes, Nn + edges]) = 2;
        edgeLengths = G.faces.areas(edges);
                            %   Contribution from each edge to
                            %   \int_\partial \Omega g_N \phi_i ds
        for e = 1:n
            bBC([nodes(2*e-1), nodes(2*e), Nn + edges(e)]) = ...
                bBC([nodes(2*e-1), nodes(2*e), Nn + edges(e)]) + ...
                edgeLengths(e).*[1/6*g(X(2*e-1,:)); 1/6*g(X(2*e,:)); 2/3*g(X(2*n + e,:))];
        end
    end
end

end