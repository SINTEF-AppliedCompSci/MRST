function [intD, intB] = faceInt(G)

           
grad_m = @(X) ...
[ones(size(X,1),1)   , zeros(size(X,1),1), zeros(size(X,1),1) ,...    %   (1,0,0)
zeros(size(X,1),1)   , ones(size(X,1),1) , zeros(size(X,1),1) ,...    %   (0,1,0)
zeros(size(X,1),1)   , zeros(size(X,1),1), ones(size(X,1),1)  ,...    %   (0,0,1)
X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1) ,...    %   (2,0,0)
X(:,2)               , X(:,1)            , zeros(size(X,1),1) ,...    %   (1,1,0)
X(:,3)               , zeros(size(X,1),1), X(:,1)             ,...    %   (1,0,1)
zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ,...    %   (0,2,0)
zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ,...    %  (0,1,1)
zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)

mx = @(X)        [X(:,1).^2/2  , ...
                     X(:,1).*X(:,2)   , ...
                     X(:,1).*X(:,3)   , ...
                     X(:,1).^3/3         , ...
                     X(:,1).^2.*X(:,2)/2 , ...
                     X(:,1).^2.*X(:,3)/2  , ...
                     X(:,1).*X(:,2).^2      , ...
                     X(:,1).*X(:,2).*X(:,3) , ...
                     X(:,1).*X(:,3).^2      ];
    
mz = @(X) ...
       [X(:,1).*X(:,3), X(:,2).*X(:,3), X(:,3).^2/2, ...
        X(:,1).^2.*X(:,3), X(:,1).*X(:,2).*X(:,3), X(:,1).*X(:,3).^2/2, ...
        X(:,2).^2.*X(:,3), X(:,2).*X(:,3).^2/2, X(:,3).^3/3];
                 
                            %   Faces and face normals.
hF = G.faces.diameters;
nF = G.faces.num;
Fc = G.faces.centroids;
faceNormals = G.faces.normals;

                            %   Edges, order by face.
edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
edges = G.faces.edges(edgeNum);
edgeSign = G.faces.edgeSign;
lengths = G.edges.lengths(edges);
normals = G.faces.edgeNormals(edgeNum,:);
nE = size(edges,1);

                            %   Nodes, ordered by face-ordered edges.
                            %   vector nodes holds all start nodes
                            %   for each edge, then all end nodes.
nodeNum = mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1);
nodes = G.edges.nodes(nodeNum);
nodes = reshape(nodes,2,[])';
nN = size(nodes,1);
nodes(edgeSign == -1,:) = nodes(edgeSign == -1,2:-1:1);
nodes = reshape(nodes,[],1);

                            %   Node coordinates and edge centroids.
                            %   Scaled for use in monomial functions.
X = G.nodes.coords(nodes,:);
Xc = G.edges.centroids(edges,:);

Xu = X;
Xcu = Xc;
X = bsxfun(@rdivide,X-repmat(rldecode(Fc, diff(G.faces.nodePos),1),2,1), ...
                      repmat(rldecode(hF, diff(G.faces.nodePos),1),2,1));
                  
Xc = bsxfun(@rdivide,Xc-rldecode(Fc, diff(G.faces.nodePos),1), ...
                        rldecode(hF, diff(G.faces.nodePos)));

%%  CALCULATE intD                                                       %%                    
                            %   Identify faces in the yz-plane, and
                            %   find corresponding edge normals, nodes
                            %   and edges centroids.
fYZ = faceNormals(:,2) == 0 & faceNormals(:,3) == 0;
eYZ = rldecode(fYZ, diff(G.faces.edgePos),1);
nYZ = repmat(rldecode(fYZ, diff(G.faces.nodePos),1),2,1);
Xx = X(~nYZ,:); nNx = size(Xx,1);
Xcx = Xc(~eYZ,:);
Xz = X(nYZ,:);  nNz = size(Xz,1);
Xcz = Xc(eYZ,:);

                            %   Scale normals by edge lengths.
normals = bsxfun(@times,normals,lengths);
xNorm = normals(~eYZ,1);
zNorm = normals(eYZ,3);

                            %   Evaluate integrals over each edge.
Ix = bsxfun(@times,(mx(Xx(1:nNx/2,:)) + mx(Xx(nNx/2+1:end,:)))/6    + ...
                    mx(Xcx)*2/3, xNorm);
Iz = bsxfun(@times,(mz(Xz(1:nNz/2,:)) + mz(Xz(nNz/2+1:end,:)))/6 + ...
                    mz(Xcz)*2/3, zNorm);
intD = zeros(nE,9);
intD(~eYZ,:) = Ix; intD(eYZ,:) = Iz;% NB: Divide by hK for grad_m!
                            
                            %   Sum for all edges of each face.
intD = mat2cell(intD, diff(G.faces.edgePos), 9);
intD = cellfun(@(X) sum(X,1), intD, 'UniformOutput', false);
intD = cell2mat(intD);
intD = bsxfun(@times, intD, hF);

%%  CALCULATE intB                                                       %%

vals = [grad_m(X); grad_m(Xc)];
vals = vals.*repmat(normals,3,9);
                            %   Fix this...
vals = [sum(vals(:,1:3),2), sum(vals(:,4:6),2),sum(vals(:,7:9),2), ...
        sum(vals(:,10:12),2),sum(vals(:,13:15),2), sum(vals(:,16:18),2), ...
        sum(vals(:,19:21),2), sum(vals(:,22:24),2), sum(vals(:,25:27),2)];
tmp = vals(G.faces.nodePos(2:end) - 1 + nN,:);
vals(nN+2:2*nN,:) = vals(nN+1:2*nN-1,:);
% valsTmp(mcolon(G.faces.nodePos(2:end-1),G.faces.nodePos(3:end)-1),:) =...
% valsTmp(mcolon(G.faces.nodePos(1:end-2),G.faces.nodePos(2:end-1)-1),:);
vals(G.faces.nodePos(1:end-1)+nN,:) = tmp;

                            %   Nodal dof integrals are found in
                            %   ... elaborate
intB = bsxfun(@rdivide,[(vals(1:nN,:) + vals(nN+1:2*nN,:))/6; vals(2*nN+1:end,:)*2/3],...
              repmat(rldecode(hF,diff(G.faces.edgePos),1),2,1));    

end