clc; clear all; close all;

P = [1/4, 0  ; ...
     3/4, 0  ; ...
     1  , 3/4; ...
     1/4, 1  ; ...
     0  , 1/2];
 


n = 8;
x = linspace(0,1,n); y= linspace(0,1,n);
[x,y] = meshgrid(x,y);
X = [x(:), y(:)];
X = X(inpolygon(X(:,1), X(:,2), P(:,1), P(:,2)),:);
np = 5;

nn = size(P,1);
a = zeros(nn,1); b = zeros(nn,1); l = ones(nn,1);
x = []; y = [];

for i = 1:nn

    p1 = P(i,:);
    p2 = P(mod(i,nn)+1,:);
    at = (p2(2)-p1(2))/(p2(1)-p1(1));
    if abs(p2(2)-p1(2)) > abs(p2(1)-p1(1))
        l(i) = 0;
        at = 1/at;
        yt = linspace(p1(2),p2(2),np*(norm(p1-p2,2)));
        bt = p2(1) - at*p2(2);
        xt = at*yt + bt;
    else
        xt = linspace(p1(1),p2(1),np*(norm(p1-p2,2)));
        bt = p2(2) - at*p2(1);
        yt = at*xt + bt;
    end
    
    a(i) = at; b(i) = bt;
    x = [x, xt(2:end)]; y = [y, yt(2:end)];

end

nn = numel(x);
for i = 1:nn
    remNodes = sum(bsxfun(@minus, X,[x(i),y(i)]).^2,2) < 1/(4*np^2);
    X = X(~remNodes,:);
end

X = [X; [x', y']];

tri = delaunay(X);

G = triangleGrid(X,tri);
G = sortEdges(G);

save('basisGrid.mat', 'G', 'P', 'a', 'b', 'l');

% P = [0, 0  ; ...
%      1, 0  ; ...
%      1, 1  ; ...
%      0, 1];
     
% tri = delaunay(P);
% nTri = size(tri,1);
% 
% G = triangleGrid(P, tri);
% G = sortEdges(G);
% G = computeGeometry(G);
% G = mrstGridWithFullMappings(G);
% 
% 
% aK = G.cells.volumes;
% n = 10;
% w = aK'/sum(aK);
% w = cumsum(w);
% 
% t = rand(n,1);
% t = diff([0,sum(bsxfun(@lt,t,w),1)]);
% 
% X = zeros(n,2);
% 
% k = 1;
% for i = 1:nTri
%     nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
%     nodes   = G.cells.nodes(nodeNum);
%     XT       = G.nodes.coords(nodes,:);
%     for j = 1:t(i)
%         at = rand(1,2);
%         xx = at(1)*(XT(2,:) - XT(1,:)) + at(2)*(XT(3,:) - XT(1,:)) + XT(1,:);
%         if ~inpolygon(xx(1), xx(2), XT(:,1), XT(:,2))
%             V = XT(2,:) + (XT(3,:)-XT(1,:));
%             xx = XT(1,:) - xx + V;
%         end
%         X(k,:) = xx;
%         k = k+1;
%     end
% end