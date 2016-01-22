function [IC, IF] = monomialCellInt(G)

%   {x, y, z, x^2, xy, xz, y^2, yz, z^}
           
int_mxx = @(X) [X(:,1).^2/2                 , ...
                X(:,1).^3/(2*3)             , ...
                X(:,1).^2.*X(:,2)/2         , ...
                X(:,1).^2.*X(:,3)/2         , ...
                X(:,1).^4/(3*4)             , ...
                X(:,1).^3.*X(:,2)/(2*3)     , ...
                X(:,1).^3.*X(:,3)/(2*3)     , ...
                X(:,1).^2.*X(:,2).^2/2      , ...
                X(:,1).^2.*X(:,2).*X(:,3)/2 , ...
                X(:,1).^2.*X(:,3).^2/2      ];

int_mxy = @(X) [X(:,1).*X(:,2), ... %  1
                X(:,1).^2.*X(:,2)/2 , ...  %   x
                X(:,1).*X(:,2).^2/2   , ...    %   y
                X(:,1).*X(:,3).*X(:,2) , ...   %   z
                X(:,1).^3.*X(:,2)/3         , ...  %   x^2
                X(:,1).^2.*X(:,2).^2/(2*2) , ...   %   xy
                X(:,1).^2.*X(:,3).*X(:,2)/2  , ... %   xz
                X(:,1).*X(:,2).^3/3      , ... %   y^2
                X(:,1).*X(:,2).^2.*X(:,3)/2 , ...
                X(:,1).*X(:,3).^2.*X(:,2)       ];

int_mx = @(X)  [X(:,1), ...
                X(:,1).^2/2  , ...
                X(:,1).*X(:,2)   , ...
                X(:,1).*X(:,3)   , ...
                X(:,1).^3/3         , ...
                X(:,1).^2.*X(:,2)/2 , ...
                X(:,1).^2.*X(:,3)/2  , ...
                X(:,1).*X(:,2).^2      , ...
                X(:,1).*X(:,2).*X(:,3) , ...
                X(:,1).*X(:,3).^2      ];
            
int_my = @(X)  [X(:,2), ...
                X(:,1).*X(:,2)  , ...
                X(:,2).^2/2   , ...
                X(:,2).*X(:,3)   , ...
                X(:,1).^2.*X(:,2)         , ...
                X(:,1).*X(:,2).^2/2 , ...
                X(:,1).*X(:,2).*X(:,3)  , ...
                X(:,2).^3/3     , ...
                X(:,2).^2.*X(:,3)/2 , ...
                X(:,2).*X(:,3).^2      ];
                        
                 
hK = G.cells.diameters;
Kc = G.cells.centroids;
nK = G.cells.num;

                            %   Faces and face normals.
faceNum = mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end)-1);
faces   = G.cells.faces(faceNum);
hF = G.faces.diameters(faces);
Fc = G.faces.centroids(faces);
faceNormals = G.faces.normals(faces,:);
nF = size(faces,1);
sign = G.faces.neighbors(faces,1) ~= rldecode((1:nK)',diff(G.cells.facePos),1);
faceNormals = bsxfun(@times,faceNormals,(-ones(nF,1)).^sign);
                            %   Unscale face normals
faceNormals = bsxfun(@rdivide, faceNormals, G.faces.areas(faces));

                            %   Edges, order by face.
edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
edges = G.faces.edges(edgeNum);
edgeSign = G.faces.edgeSign(edgeNum);
lengths = G.edges.lengths(edges);
edgeNormals = G.faces.edgeNormals(edgeNum,:);
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

                            %   Node coordinates and Gauss Lobatto
                            %   quadrature points:
                            %   X(1:nN,:)-Xq(1:nN,:)-Xq(nN+1:2*nN,:)-X(nN+1:2*nN,:)
X = G.nodes.coords(nodes,:);
edgeVec = X(nN+1:2*nN,:) - X(1:nN,:);
Xq = [X(1:nN,:) + edgeVec*(0.5 - sqrt(1/20)); ...
      X(1:nN,:) + edgeVec*(0.5 + sqrt(1/20))]; ...
                            
                            %   Find number of face-nodes per cell
diffVec = diff(G.faces.nodePos);
cellFaceNodes = cellfun(@(X) sum(X,1), ...
                mat2cell(diffVec(faces), diff(G.cells.facePos), 1));
  
Xu = X;

                            %   Scale coordinates for use in monomials.
X = bsxfun(@rdivide,X-repmat(rldecode(Kc,cellFaceNodes,1),2,1), ...
                      repmat(rldecode(hK, cellFaceNodes,1),2,1));
Xq = bsxfun(@rdivide,Xq-repmat(rldecode(Kc,cellFaceNodes,1),2,1), ...
                        repmat(rldecode(hK, cellFaceNodes,1),2,1));
                    
%%  CALCULATE intD                                                       %%

                           %   Identify faces in the yz-plane, and
                            %   find corresponding edge normals, nodes
                            %   and edges centroids.
diffVec = diff(G.faces.edgePos);
fYZ = faceNormals(:,2) == 0 & faceNormals(:,3) == 0;
eYZ = rldecode(fYZ, diffVec(faces),1);
diffVec = diff(G.faces.nodePos);
nYZ = repmat(rldecode(fYZ, diffVec(faces),1),2,1);
Xx = X(~nYZ,:); nNx = size(Xx,1);
Xqx = Xq(~nYZ,:);
Xy = X(nYZ,:);  nNy = size(Xy,1);
Xqy = Xq(nYZ,:);
                            %   Scale normals by edge lengths.
edgeNormals = bsxfun(@times,edgeNormals,lengths);
                            %   Product of first component of face normal
                            %   and edge normal.

diffVec = diff(G.faces.edgePos);
normProdx = edgeNormals(~eYZ,1).*rldecode(faceNormals(~fYZ,1), diffVec(faces(~fYZ)),1);
normPrody = edgeNormals(eYZ,2).*rldecode(faceNormals(fYZ,1), diffVec(faces(fYZ)),1);

                            %   Evaluate integrals over each edge.
intX = bsxfun(@times, ...
              (int_mxx(Xx(1:nNx/2,:))  + int_mxx(Xx(nNx/2+1:nNx,:)))/12 + ...
              (int_mxx(Xqx(1:nNx/2,:)) + int_mxx(Xqx(nNx/2+1:nNx,:)))*5/12, ...
               normProdx);
intY = bsxfun(@times,(int_mxy(Xy(1:nNy/2,:))  + int_mxy(Xy(nNy/2+1:nNy,:)))/12 + ...
                     (int_mxy(Xqy(1:nNy/2,:)) + int_mxy(Xqy(nNy/2+1:nNy,:)))*5/12, ...
                     normPrody);
                            %   Assemble all integrals
IC = zeros(nE,10); IC(~eYZ,:) = intX; IC(eYZ,:) = intY;
                            %   Sum for all edges of each face.

IC = mat2cell(IC, cellFaceNodes, 10);
IC = cellfun(@(X) sum(X,1), IC, 'UniformOutput', false);
IC = bsxfun(@times, cell2mat(IC), hK.^2);

faceNodes = diff(G.faces.nodePos);

Xq = G.edges.centroids(edges,:);
Xq = bsxfun(@rdivide,Xq-rldecode(Kc,cellFaceNodes,1), ...
                        rldecode(hK,cellFaceNodes,1));

Xqx = Xq(~eYZ,:);
Xqy = Xq(eYZ,:);

                            %   Evaluate integrals over each edge.
intX = bsxfun(@times, ...
              (int_mx(Xx(1:nNx/2,:))  + int_mx(Xx(nNx/2+1:nNx,:)))/6 + ...
               int_mx(Xqx)*2/3, ...
               edgeNormals(~eYZ,1));
intY = bsxfun(@times,(int_my(Xy(1:nNy/2,:))  + int_my(Xy(nNy/2+1:nNy,:)))/6 + ...
                      int_my(Xqy(1:nNy/2,:))*2/3, ...
                     edgeNormals(eYZ,2));

                            %   Assemble all integrals
IF = zeros(nE,10); IF(~eYZ,:) = intX; IF(eYZ,:) = intY;
                 
IF = mat2cell(IF, faceNodes(faces), 10);
IF = cellfun(@(X) sum(X,1), IF, 'UniformOutput', false);
IF = bsxfun(@times, cell2mat(IF), hK);

%   Speed improvements: First columns of IC and IF are just volume and
%   area.

end