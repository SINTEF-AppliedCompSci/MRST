function G = voronoi2mrstGrid3D(V, C)
% Transform Voronoi diagram stored in a Qhull grid structure to a MRST grid
% structure.
%
% SYNOPSIS:
%   G = voronoi2mrstGrid3D(V,C)
%
% PARAMETERS:
%   V         A nx3 array containing the vertices of the Voronoi diagram,
%             as obtained from [V, C] = voronoin(pts)
%   C         A cell array where each element coresponds to one cell,
%             as obtained from [V, C] = voronoin(pts)
%
% RETURNS:
%   G                - Valid MRST grid definition.  
%
% EXAMPLE:
% dt = 0.2;
% [X,Y,Z] = ndgrid(0:dt:1);
% X(1:2:end) = X(1:2:end) + dt/2;
% Y(1:2:end) = Y(1:2:end) + dt/2;
% pts = [X(:), Y(:),Z(:)];
% [V,C] = voronoin(pts);
% G = voronoi2mrstGrid3D(V,C);
% plotGrid(G);
% axis equal
%
% SEE ALSO:
%   compositePebiGrid2D, pebi, clippedPebi3D, clipGrid.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

assert(size(V, 2) == 3, ...
      ['Function ''%s'' is only supported in three ', ...
       'space dimensions.'], mfilename);
     
if isinf(V(1,1))
  warning(['Found cells extending to infinity. These cells will be ', ...
           'removed from the converted grid'])
  rem = cellfun(@(c) any(isinf(V(c,1))), C);
  C   = C(~rem);
  V   = V(2:end,:);
  C   = cellfun(@(c) c-1, C,'un',false);
end

% Create mapping from cells to nodes.
cell2Node    = cumsum([1; cellfun(@numel, C)])';
activeVertex = cell2mat(C');
[activeVertex, ~, C] = unique(activeVertex);
V            = V(activeVertex,:);

% Set number of cells
G.cells.num  = numel(cell2Node)-1;

% Find half faces for all cells
facePos    = ones(G.cells.num+1,1);
hf         = [];      
hf2NodePos = [1]; 

for i = 1:G.cells.num
    % Calculate convex hull
    H = convhull(V(C(cell2Node(i):cell2Node(i+1)-1),:),'simplify',true);
    hull         = C(cell2Node(i)-1+H);
    % Merge triangle faces into polygons
    [hull, localPos] = remParFaces(V, hull);
    hf           = [hf; hull];
    localPos = unique(localPos);
    hf2NodePos   = [hf2NodePos; hf2NodePos(end)-1 + localPos(2:end)];
    facePos(i+1) = numel(hf2NodePos-1);

end
G.cells.facePos = facePos;

% Find faces from the half faces.
[nodes, nodePos, ic] = uniqueFace(hf,hf2NodePos);

G.faces.nodePos = nodePos;
G.faces.nodes   = reshape(nodes', [], 1);
G.cells.faces   = ic;
G.nodes.coords  = V;
G.nodes.num     = size(V,1);  
G.faces.num     = max(G.cells.faces);

% Set neighbors
cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
G.faces.neighbors = zeros(G.faces.num,2);
for i = 1:G.faces.num
    neigh = G.cells.faces==i;
    if sum(neigh)==2
        G.faces.neighbors(i,[1,2]) = cellNo(neigh);
    else
        G.faces.neighbors(i,:) = [cellNo(neigh),0];
    end
end

% Set grid info
G.type    = { mfilename };
G.griddim = 3;
end


function [faceNodes,nodePos, ic] = uniqueFace(hfNodes, hf2node)
% Finds the faces given half faces.

% Find all half faces with equal number of nodes
faceSize  = diff(hf2node);
[~,ias,~] = unique(faceSize);

faceNodes = [];
nodePos   = [1];
ic = zeros(size(hf2node,1)-1,1);
for i = 1:numel(ias) % For all half faces with equal number of nodes
    % Find the indexes of the half face nodes
    testPos  = faceSize(ias(i))==faceSize;
    fromFace = hf2node([testPos;false]);
    toFace   = hf2node([false;testPos]) - 1;
    nodeID   = arrayfun(@(l,r) (l:r), fromFace, toFace, 'un', 0)';
    nodeID   = cell2mat(nodeID);
    tempFace = reshape(hfNodes(nodeID),faceSize(ias(i)),[]);
    tempFace = tempFace';
    if isempty(tempFace)
      continue
    end
    % Half faces with the same nodes are one face
    [~,ia2,ic2] = unique(sort(tempFace,2),'rows');
    temp        = tempFace(ia2,:)';
    faceNodes   = [faceNodes;temp(:)];
    ic(testPos) = ic2+numel(nodePos)-1;
    nodePos     = [nodePos; ...
                   nodePos(end)+cumsum(repmat(faceSize(ias(i)),[size(temp,2),1]))];
end
end

function [newHull, nodePos] = remParFaces(V, hull)
% Merge all faces in convex hull that have equal normals

newHull = [];
nodePos = [1];
% Calculate normals
n       = calcNormals(hull, V);

% This might be necesary
%    e = any(isnan(n),2);
%    n = n(~e,:);
%    hull = hull(~e,:);

while ~isempty(hull)
    % Find face normals that are equal to first face normal in stack
    parFace = n*(n(1,:)')> 1 - 50*eps;
    % Merge parallel faces
    tmp     = mergeFaces(V, hull, parFace, n(1,:)');
    nf      = numel(tmp);
    if nf ==0
      hull    = hull(~parFace,:);
      n       = n(~parFace,:);
      continue
    end
    newHull = [newHull; tmp];
    nodePos = [nodePos; nodePos(end) + nf];
    % Update hull
    hull    = hull(~parFace,:);
    n       = n(~parFace,:);
end
end


function [merged] = mergeFaces(V, H, F, n)
% Merge faces nodes in counterclockwise direction

% Find unique node index
merged     = unique(reshape(H(F,:),[],1));
% shift coordinate system
id         = find(F);
x0         = mean(V(H(id(1),:),:));
VC         = bsxfun(@minus, V, x0);
% Create new basis
basis      = VC(H(id(1),1:2),:)';
basis(:,1) = basis(:,1)/norm(basis(:,1),2);
basis(:,2) = basis(:,2)-(basis(:,1)'*basis(:,2))*basis(:,1);
basis(:,2) = basis(:,2)/norm(basis(:,2),2);

% find coordinates in new basis
VB         = (basis\VC(merged,:)')';

% Sort nodes based on the angle
theta      = atan2(VB(:,2),VB(:,1));
[~,i]      = sort(theta);
merged     = merged(i);

% Remove colinear points
merged = [merged; merged(1:2)];
j = 1;
k = 2;
rem = false(size(merged,1),1);

while j<size(merged,1)-1
    if isColinear(V([merged(j);merged(k:k+1)],:));
       rem(k) = true;
       k = k+1;
    else
        j = k;
        k = j+1;
    end
    if  k==size(merged,1)
        break
    end
end
merged = merged([~rem(end-1);~rem(2:end-2);false;false],:);
if numel(merged)<=2 %face is a line
  merged = [];
end
end
