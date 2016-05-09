function [A,b] = VEM3D_bc(G, A, b, bc, k)
%--------------------------------------------------------------------------
%   Incorporates boundary conditions in stiffness matrix A and load term b,
%   obtained using a kth order VEM.
%
%   SYNOPSIS:
%       [A, b] = VEM2D_bc(G, A, b, bc, k)
%
%   DESCRIPTION:
%       Incorporates boundary conditions in stiffness matrix A and load
%       term b, obtained using a kth order VEM. Dirichlet conditions are
%       set directly in the load term, and the corresponding row of A is
%       changed to the corresponding row of the identity matrix. Neumann
%       conditions adds a boundary integral to the load term, which is
%       evaluated using a k+1 accurate quadrature rule. See [1] for
%       details.
%
%   REQUIRED PARAMETERS:
%       G   - MRST grid.
%       A   - Global VEM stiffness matrix.
%       b   - Global VEM load term.
%       bc  - Boundary condition struct constructed using VEM2D_addBC.
%       k   - Method order.
%
%   RETURNS:
%       A   - Global stiffness matrix with boundary condition
%             modifications.
%       b   - Global load term, with boundary condition modifications.    
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

assert(k == 1 | k == 2); 

faces = bc.faces(strcmp(bc.type,'flux'));
func = bc.func(strcmp(bc.type,'flux'));
nF = numel(faces);
if ~isempty(faces)
    [m, ~, ~] = retrieveMonomials(2,k);
    [Xq, w, ~, vol] = triangleQuadRule(k+1);
    nq = size(Xq,1);
end
    
for i = 1:nF
    
    F = faces(i);
    g = func{i};

    %   Edge data for face F.

    edgeNum     = G.faces.edgePos(F):G.faces.edgePos(F+1)-1;
    edges       = G.faces.edges(edgeNum);
    
    %   Node data for each edge of each face.
    
    nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
    nodes   = G.edges.nodes(nodeNum);
    nodes   = reshape(nodes,2,[])';
    nN      = size(nodes,1);
    nodes(G.faces.edgeSign(edgeNum) == -1,:) ...
            = nodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
    nodes   = nodes(:,1);
    
    if k == 1
        X = G.nodes.coords(nodes,:);
        dofVec = nodes';
    elseif k == 2
        X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:)];
        dofVec = [nodes', edges' + G.nodes.num, ...
                                            F + G.nodes.num + G.edges.num];
    end
    xx = X;
    
    T = G.faces.localCoords(G.faces.TPos(F):G.faces.TPos(F+1)-1,:);
    X = (bsxfun(@minus, X, G.faces.centroids(F,:)))*T;

    tri = delaunay(X(1:nN,:));
    nTri = size(tri,1);
    
    %   Construct map from refrence triangle to triangles in triangulation.
    
    bA = X(tri(:,1),:);
    M = X(tri(:,2:end),:) - repmat(bA,2,1);
    M = M(mcolon(1:nTri,2*nTri,nTri),:);
    M = mat2cell(M,2*ones(nTri,1),2);
    detA = cellfun(@(X) abs(det(X)), M);

    %   Map quadrature points from reference triangle to triangels.
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, M, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);
                            
    %   Evaluate functions at quadrature points.
    
    Xmon = Xhat/G.faces.diameters(F);

    PNFstar = G.faces.PNstarT(G.faces.PNstarPos(F):...
                               G.faces.PNstarPos(F+1)-1,:)';
    Xhat = bsxfun(@plus, Xhat*T', G.faces.centroids(F,:));
    
%     plot3(xx(:,1), xx(:,2),xx(:,3), '.');
%     hold on
%     plot3(Xhat(:,1), Xhat(:,2), Xhat(:,3), 'o');
%     
    vals = bsxfun(@times,m(Xmon)*PNFstar,g(Xhat));
    
    %   Multilply by wheights and determinants.
    
    detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
    vals = bsxfun(@times,vals,detAw);

    b(dofVec) = b(dofVec) + sum(vals,1)';
    
end

faces = bc.faces(strcmp(bc.type, 'pressure'));
func = bc.func(strcmp(bc.type,'pressure'));
nF = numel(faces);

for i = 1:nF
    
    F = faces(i);
    g = func{i};
    
    nodeNum = G.faces.nodePos(F):G.faces.nodePos(F+1)-1;
    nodes   = G.faces.nodes(nodeNum);

    if k == 1
        X = G.nodes.coords(nodes,:);
        gChi = g(X);
        dofVec = nodes';
    elseif k == 2
        edgeNum = G.faces.edgePos(F):G.faces.edgePos(F+1)-1;
        edges   = G.faces.edges(edgeNum);
        int = polygonInt3D(G, F, g, k+1)./G.faces.areas(F);
        X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:)];
        gChi = [g(X); int];
        dofVec = [nodes', edges' + G.nodes.num, ...
                                            F + G.nodes.num + G.edges.num];
    end
    b(dofVec) = gChi;
end

nodeNum = mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1);
nodes   = G.faces.nodes(nodeNum);
if k == 1
    dofVec = nodes';
elseif k == 2
    edgeNum = mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1);
    edges = G.faces.edges(edgeNum);
    dofVec = [nodes', edges' + G.nodes.num, faces' + G.nodes.num + G.edges.num];
end

N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2 ...
                                                 + G.cells.num*k*(k^2-1)/6;
I = spdiags(ones(N,1),0,N,N);
A(dofVec,:) = I(dofVec,:);