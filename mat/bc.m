% delta = 0.001;
% gD = @(X) 10*ones(size(X,1),1);
% gN = @(X) zeros(size(X,1),1);
% f = @(X) 1/sqrt(2*delta*pi).*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2)./(2*delta));
% 
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bNeu = zeros(numel(boundaryEdges),1);
% for e = 1:numel(boundaryEdges)
%     X = G.nodes.coords(G.faces.nodes(G.faces.nodePos(e):G.faces.nodePos(e +1)-1));
%     if X(:,1) < 0.001
%         bNeu(e) = 1;
%     end
% end
% bc = struct('bcFunc', {{gN,gD}}, 'bcFaces', {{boundaryEdges(bNeu == 1), boundaryEdges(bNeu == 0)}}, 'bcType', {{'neu', 'dir'}});

% gD = @(X) 1/4.*((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2);
% f = @(X) -ones(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});


% g = @(X) 1/6.*((X(:,1)-0.5).^3 + (X(:,2)-0.5).^3);
% f = @(X) -(X(:,1)-0.5) - (X(:,2)-0.5);

% g = @(X) X(:,2).*(1-X(:,2)).*X(:,1).^3;
% f = @(X) -6.*X(:,1).*X(:,2).*(1-X(:,2)) + 2.* X(:,1).^3;

% gD = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
% f = @(X) zeros(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

gD = @(X) 1/2.*(X(:,2).^2 - X(:,1).^2);
gN = @(X) X(:,2);
f = @(X) zeros(size(X,1),1);
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bNeu = zeros(numel(boundaryEdges),1);
bDir = zeros(numel(boundaryEdges),1);
for e = 1:numel(boundaryEdges)
    X = G.faces.centroids(boundaryEdges(e),:);
    if X(:,1) < 0 + eps 
        bDir(e) = 1;
    end
    if X(:,1) > 1 - eps
        bDir(e) = 1;
    end
    if X(:,2) < 0 + eps 
        bNeu(e) = 1;
    end
    if X(:,2) > 1-eps
        bNeu(e) = 1;
    end
    
end

bc = struct('bcFunc', {{gN,gD}}, 'bcFaces', {{boundaryEdges(bNeu == 1), boundaryEdges(bDir == 1)}}, 'bcType', {{'neu', 'dir'}});
