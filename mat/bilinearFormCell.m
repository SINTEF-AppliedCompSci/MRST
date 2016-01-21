function I = bilinearFormCell(G)

grad_m = @(X) ...
[ones(size(X,1),1)   , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (1,0,0)
zeros(size(X,1),1)   , ones(size(X,1),1) , zeros(size(X,1),1) ;...    %   (0,1,0)
zeros(size(X,1),1)   , zeros(size(X,1),1), ones(size(X,1),1)  ;...    %   (0,0,1)
X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (2,0,0)
X(:,2)               , X(:,1)            , zeros(size(X,1),1) ;...    %   (1,1,0)
X(:,3)               , zeros(size(X,1),1), X(:,1)             ;...    %   (1,0,1)
zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ;...    %   (0,2,0)
zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ;...    %  (0,1,1)
zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)


edgeNum = G.faces.edgePos(1):G.faces.edgePos(2)-1;
faceEdges = G.faces.edges(edgeNum);
edgeNormals = G.faces.edgeNormals(edgeNum,:);
    
nodeNum = mcolon(G.edges.nodePos(faceEdges),G.edges.nodePos(faceEdges+ 1)-1);
faceNodes = G.edges.nodes(nodeNum);
faceNodes = reshape(faceNodes,2,[])';
nNF = size(faceNodes,1);


NF = 2*nNF+1;
faceNodes(G.faces.edgeSign(edgeNum) == -1,:) = ...
    faceNodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
faceNodes = reshape(faceNodes,[],1);
vec1 = (G.nodes.coords(faceNodes(nNF+1),:) - G.nodes.coords(faceNodes(1),:))';
                        %   NB: faceNormal must be in correct
                        %   direction, as done by faceData.
vec2 = cross(G.faces.normals(1,:)',vec1);
T = [vec1,vec2];

X = G.nodes.coords(faceNodes,:);
n = G.faces.normals(1,:)*T;

I = grad_m(X)*T;

    

end 