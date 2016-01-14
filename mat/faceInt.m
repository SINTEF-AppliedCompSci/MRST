function [ID, IB] = faceInt(G)

    
m =      @(X) [X(:,1)               , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,3)               , ...   %   (0,0,1)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,1).*X(:,3) , ...   %   (1,0,1)
               X(:,2).^2 , ...   %   (0,2,0) 
               X(:,2).*X(:,3) , ...   %   (0,1,1)
               X(:,3).^2          ];      %   (0,0,2)
grad_m = @(X) [ones(size(X,1),1)    , zeros(size(X,1),1),  zeros(size(X,1),1);...    %   (1,0,0)
               zeros(size(X,1),1)   , ones(size(X,1),1) ,  zeros(size(X,1),1);...    %   (0,1,0)
               zeros(size(X,1),1)   , zeros(size(X,1),1),  ones(size(X,1),1);...    %   (0,0,1)
               X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1);...    %   (2,0,0)
               X(:,2)               , X(:,1)            , zeros(size(X,1),1);...    %   (1,1,0)
               X(:,3)               , zeros(size(X,1),1), X(:,1)   ;...    %   (1,0,1)
               zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1);...    %   (0,2,0)
               zeros(size(X,1),1)   , X(:,3)            , X(:,2)   ;...    %  (0,1,1)
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
                 
    
    hF = G.faces.diameters;
    nF = G.faces.num;
    Fc = G.faces.centroids;
    faceNormals = G.faces.normals;
    
    edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
    edges = G.faces.edges(edgeNum);
    lengths = G.edges.lengths(edges);
    normals = G.faces.edgeNormals(edgeNum,:);
    
    nodeNum = mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1);
    nodes = G.edges.nodes(nodeNum);
    
    X = [G.nodes.coords(nodes,:);
         G.edges.centroids(edges,:)];
    C = rldecode([Fc; Fc], ...
         [2.*diff(G.faces.nodePos);diff(G.faces.nodePos)]);
    hFNodes = rldecode([hF;hF], ...
         [2.*diff(G.faces.nodePos);diff(G.faces.nodePos)],1);
    
     
    
    X = bsxfun(@rdivide,(X - C),hFNodes);
    clear hFNodes C
    nV = size(nodes,1);
    
    fInYZ = faceNormals(:,2) == 0 & faceNormals(:,3) == 0;
    eInYZ = rldecode(fInYZ,...
            diff(G.faces.edgePos),1);
    vInYZ = rldecode([fInYZ; fInYZ], ...
         [2.*diff(G.faces.nodePos);diff(G.faces.nodePos)],1);
    Xx = bsxfun(@times,X,~vInYZ);
    Xz = bsxfun(@times,X,vInYZ);
    
    Xx = X(1:nV,:);
    Xx = Xx(~vInYZ(1:nV),:);
    Xz = X(1:nV,:);
    Xz = Xz(vInYZ(1:nV),:);
    XxC = X(nV+1:end,:);
m =      @(X) [X(:,1)               , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,3)               , ...   %   (0,0,1)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,1).*X(:,3) , ...   %   (1,0,1)
               X(:,2).^2 , ...   %   (0,2,0) 
               X(:,2).*X(:,3) , ...   %   (0,1,1)
               X(:,3).^2          ];      %   (0,0,2)
    XxC = XxC(~vInYZ(nV+1:end),:);
    XzC = X(nV+1:end,:);
    XzC = XzC(vInYZ(nV+1:end),:);
    
    aN = bsxfun(@times,normals,lengths);
    aNx = aN(~eInYZ,1);
    aNz = aN(eInYZ,3);
    
    Ix = bsxfun(@times,(mx(Xx(1:2:end-1,:)) + mx(Xx(2:2:end,:)))/6    + ...
                        mx(XxC)*2/3, aNx);
    
    Iz = bsxfun(@times,(mz(Xz(1:2:end-1,:)) + mz(Xz(2:2:end,:)))/6 + ...
                         mz(XzC)*2/3, aNz);
    nE = size(edges,1);
                     
    ID = zeros(nE,9);
    ID(~eInYZ,:) = Ix;
    ID(eInYZ,:) = Iz;
    
    
    diffVec = diff(G.faces.edgePos);
    it = 1;
    for i = 1:nF
        ID(i,:) = sum(ID(it:it + diffVec(i)-1,:),1);
        it = it + diffVec(i);
    end
    ID = bsxfun(@times, ID(1:nF,:), hF);
                            % NB: Divide by hK for grad_m!
    IB = sum([grad_m(X(1:2:nV-1,:))/6; grad_m(X(2:2:nV,:))/6; ...
               grad_m(X(nV+1:end,:))*2/3].*repmat(aN,9,1),2);
    IB = reshape(IB,nE,[]);
    it = 1;
    for i = 1:nF
        IB(i,:) = sum(IB(it:it + diffVec(i)-1,:),1);
        it = it + diffVec(i);
    end
    IB
    
end