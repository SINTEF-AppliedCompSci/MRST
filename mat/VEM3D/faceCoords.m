function G = faceCoords(G)

nF          = G.faces.num;
faceNormals = G.faces.normals;

% edgeNum     = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
% edges       = G.faces.edges(edgeNum);
% edgeNormals = G.faces.edgeNormals(edgeNum,:);

nodes = G.faces.nodes(mcolon(G.faces.nodePos(1:end-1),...
                                              G.faces.nodePos(2:end)-1));             
XN = G.nodes.coords(nodes,:);
numFaceNodes = diff(G.faces.nodePos);
nN = numel(nodes);

% X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:);];

%%  MAP FROM GLOBAL TO LOCAL COORDINATES                                 %%

%   Build local coordinate systems. x -> x*T' + b.



vec1 = (XN(G.faces.nodePos(1:end-1)+1,:) - XN(G.faces.nodePos(1:end-1),:));
vec1 = bsxfun(@rdivide, vec1, sqrt(sum(vec1.^2,2)));
vec2 = cross(faceNormals,vec1,2);
vec2 = bsxfun(@rdivide, vec2, sqrt(sum(vec2.^2,2)));
vec1 = vec1'; vec2 = vec2';

C = sparseBlockDiag([vec1(:),vec2(:)],3*ones(1,nF),1);
d    = G.faces.centroids;

XN = XN-rldecode(d, numFaceNodes);
XN = sparseBlockDiag(XN, numFaceNodes,1);
XFN   = XN*C;

edgeNum = mcolon(G.faces.edgePos(1:end-1), G.faces.edgePos(2:end)-1);
edges = G.faces.edges(edgeNum);
numFaceEdges = diff(G.faces.edgePos);
XE = G.edges.centroids(edges,:);
nE = numel(edges);
XE = XE-rldecode(d, numFaceNodes);
XE = sparseBlockDiag(XE, numFaceEdges,1);
XFE   = XE*C;

NE = G.faces.edgeNormals(edgeNum,:);
NE = sparseBlockDiag(NE, numFaceEdges,1);
NFE   = NE*C;

XFN = squeezeBlockDiag(XFN, numFaceNodes, nN, 2);
XFE = squeezeBlockDiag(XFE, numFaceEdges, nE, 2);
NFE = squeezeBlockDiag(NFE, numFaceEdges, nE, 2);
C = squeezeBlockDiag(C, numFaceNodes, 3*nF, 2);
CPos = 1:3:3*nF;

G.faces.faceCoords = XFN;
G.faces.edgeCentroids = XFE;
G.faces.edgeNormalsC = NFE;
G.faces.C = C;
G.faces.CPos = CPos;