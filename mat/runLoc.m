clc; clear all; close all;
run('../../matlab/project-mechanics-fractures/mystartup.m')

                             %   Build cartesian grid and compute geometry.
nx = 3; ny = 3  ;            %   Grid dimensions.
G = cartGrid([nx, ny]);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);
Nc = G.cells.num;            %   Number of cells in G
Ne = G.faces.num;
Nn = G.nodes.num;

c = 5;
f = 1;
                            %   Get coordinates of vertices of cell c.
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
                            %   Find edge midpoints.
    faceNum = G.cells.facePos(c) : G.cells.facePos(c+1)-1;
    faces = G.cells.faces(faceNum);
    Xmid = G.faces.centroids(faces,:);
                            %   Find boundary edges.
    neighbors = G.faces.neighbors(faces,:);
    boundaryEdges = sum(neighbors == 0, 2);
                            %   Find edge normals, fix orientation.
    normals = G.faces.normals(faces,:);
    m = (-ones(length(normals),1)).^(neighbors(:,1) ~= c);
    normals = [m,m].*normals;
    edgeLengths = G.faces.areas;
                            %   Find volume.
    vol = G.cells.volumes(c);
                            %   Calculate local stiffness matrix.
    X = X([2:end,1],:);
    Xmid = Xmid([2:end,1],:);
    normals = normals([2:end,1],:);
    edgeLengths = edgeLengths([2:end,1]);
    
    [Sl, bl] = locS(X, Xmid, edgeLengths, normals, boundaryEdges, vol, f,sqrt(2));
    
    hK = 5;
    X = [0,0;3,0;3,2;3/2,4;0,4]
    n = size(X,1);
    Xmid = (X(1:n, :) + X([2:end, 1], :))/2
    [~, XB] = baric(X);
    vol = 21/2;
    boundaryEdges = zeros(1,n);
    normals = zeros(n,2);
    for i = 1:n
         nVec = [ X(mod(i,n)+1,2) - X(i,2) , - X(mod(i,n)+1,1) + X(i,1) ];
         nVec = nVec./norm(nVec);
         normals(i,:) = nVec;
    end
    edgeLengths = zeros(n,1);
    for i = 1:n
        edgeLengths(i) = sqrt((-X(mod(i,n)+1,1) + X(i,1))^2 + (X(mod(i,n)+1,2) - X(i,2))^2);
    end
    
    [Sl, bl] = locS(X, Xmid, edgeLengths, normals, boundaryEdges, vol, f, 5);
    

% Ndof = Nc + Nn + Ne;
% S = sparse(Ndof,Ndof);
% f = @(X) sin(X(:,1).*pi).*sin(X(:,2).*pi);
% b = zeros(Ndof,1);

% for c = 1:Nc
% 
%     %  CALUCLATE LOCAL STIFFNESS MATRIX FOR CELL c                      %%
%                               Get coordinates of vertices of cell c.
%     nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
%     nodes = G.cells.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:);
%                               Find edge midpoints.
%     faceNum = G.cells.facePos(c) : G.cells.facePos(c+1)-1;
%     faces = G.cells.faces(faceNum);
%     Xmid = G.faces.centroids(faces,:);
%                               Find boundary edges.
%     neighbors = G.faces.neighbors(faces,:);
%     boundaryEdges = sum(neighbors == 0, 2);
%                               Find edge normals, fix orientation.
%     normals = G.faces.normals(faces,:);
%     m = (-ones(length(normals),1)).^(neighbors(:,1) ~= c);
%     normals = [m,m].*normals;
%                               Find volume.
%     vol = G.cells.volumes(c);
%                               Calculate local stiffness matrix.
%     [Sl, bl] = locS(X, Xmid, normals, boundaryEdges, vol, f);
% 
%     
%     %  MAP LOCAL STIFFNESS MATRIX TO GLOBAL STIFFNESS MATRIX            %%
% 
%     dofVec = [nodes', faces + Nn, Nn + Ne + c];
%     
%     S(dofVec, dofVec) = S(dofVec, dofVec) + Sl;
%     b(dofVec) = b(dofVec) + bl;
% end
% 
% sol = S\b;
% X = [G.nodes.coords; G.faces.centroids; G.cells.centroids];
% [x,y] = meshgrid(0:0.5:nx, 0:0.5:ny);
% z = griddata(X(:,1), X(:,2), sol, x, y);
% surf(x,y,z)
% 
%    
%     dofvec=mcolon(2*(nodes-1)+1,2*(nodes-1)+2)
%     iindl=repmat(dofvec',1,numel(dofvec));
%     jindl=repmat(dofvec,numel(dofvec),1);
%     nnzl = numel(iindl);
%     iglob(pos+1:pos+nnzl) = iindl;
%     jglob(pos+1:pos+nnzl) = jindl;
%     mglob(pos+1:pos+nnzl) = Sl;
%     pos = pos+nnzl;