function [bcDof, bBC] = VEM3D_bc(G, BC, k)
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
nE = G.edges.num;           %   Number of edges.
nF = G.faces.num;
nK = G.cells.num;           %   Number of cells.

N = nN + nE*(k-1) + nF*k*(k-1)/2 + nK*k*(k^2-1)/6;        %   Number of dofs.

bcFunc = BC.bcFunc;         %   Boundary condition functions.
bcFaces = BC.bcFaces;       %   Vectors of boundary edges.
bcType = BC.bcType;         %   Boundary condition types
nBC = numel(bcFunc);        %   Number of boundary conditions.
bcDof = zeros(N, 1);     %   Dof map for boundary conditions.
bBC = zeros(N, 1);       %   b vector for boundary conditions.
for b = 1:nBC
    g = bcFunc{b};          %   BC funciton.
    faces = bcFaces{b};     %   BC faces.
    type = bcType{b};       %   BC type.
    
    nodeNum = mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1);
    nodes = G.faces.nodes(nodeNum);
    
    if k == 1
        
        X = G.nodes.coords(nodes,:);
        gChi = g(X);
        dofVec = nodes';
        
    
    elseif k == 2
    
        edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
        edges = G.faces.edges(edgeNum);
    
        I = polygonInt3D(G, faces, g, 3)./G.faces.areas(faces);
    
        X = [G.nodes.coords(nodes,:)    ; ...
             G.edges.centroids(edges,:)];
         
        gChi = [g(X); I];
         
        dofVec = [nodes', edges' + nN, faces' + nN + nE*(k-1)];
    end
     
    switch type
        
        case 'dir'
            bcDof(dofVec) = 1;
            bBC(dofVec) = gChi;
        case 'neu'
%             bcDof([nodes, Nn + edges]) = 2;
%             edgeLengths = G.faces.areas(edges);
    end
            
%             
%     if strcmp(type, 'dir')
%                             %   Apply Dirichelt BC's.
%         
%     elseif strcmp(type, 'neu')
%                             %   Apply Neuman BC's.
%         ;
%                             %   Contribution from each edge to
%                             %   \int_\partial \Omega g_N \phi_i ds
%         for e = 1:n
%             bBC([nodes(2*e-1), nodes(2*e), Nn + edges(e)]) = ...
%                 bBC([nodes(2*e-1), nodes(2*e), Nn + edges(e)]) + ...
%                 edgeLengths(e).*[1/6*g(X(2*e-1,:)); 1/6*g(X(2*e,:)); 2/3*g(X(2*n + e,:))];
%         end
%     end
end

end