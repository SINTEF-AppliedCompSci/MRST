function I = faceInt(G)

    m_int = @(X,hK)  hK.*[X(:,1).^2/2  , ...
                     X(:,1).*X(:,2)   , ...
                     X(:,1).*X(:,3)   , ...
                     X(:,1).^3/3         , ...
                     X(:,1).^2.*X(:,2)/2 , ...
                     X(:,1).^2*X(:,3)/2  , ...
                     X(:,1).*X(:,2).^2      , ...
                     X(:,1).*X(:,2).*X(:,3) , ...
                     X(:,1).*X(:,3).^2/2      ];
    
    nV = G.nodes.num;
    hK = G.faces.diameters;
    Fc = G.faces.centroids;
    
    edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
    edges = G.faces.edges(edgeNum);
    
    nodeNum = mcolon(G.faces.nodePos(1:end-1), G.faces.nodePos(2:end)-1);
    nodes = G.faces.nodes(nodeNum);
    
    X = [G.nodes.coords(nodes,:);
         G.edges.centroids(edges,:)];
    X = (X - rldecode([Fc; Fc], repmat(diff(G.faces.nodePos),2,1),1))./...
        repmat(rldecode([hK;hK],repmat(diff(G.faces.nodePos),2,1),1),1,3);
    
%     figure;
%     axis([-1,2,-1,2,-1,2]);
%     view(3);
%     hold on
%     for i = 1:size(X,1)
%         plot3(X(i,1), X(i,2), X(i,3),'o')
%         pause
%     end
%     
%     G.faces.faceDofs
    
    
end