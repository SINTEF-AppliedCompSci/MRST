function [V, C, symV] = clipGrid(dt,bound,varargin)
% Calculate the intersection of a Voronoi diagram and a surface in 3D
%
% SYNOPSIS:
%   [V, C, symV] = clipGrid(dt, bound)
%   [...] = clipGrid(...,'Name1',Value1,'Name2',Value2)
% PARAMETERS:
%   dt       A Delaunay triangulation of the Voronoi sites. dt must contain
%            the following parameters, see the matlab function 
%            delaunayTriangulation for further information:
%               - dt.Points
%               - dt.Edges
%
%   bound    A triangulation of the surface. bound must contain the
%            following information, as described in the Matlab function 
%            delaunayTriangulation:
%               - bound.Points
%               - bound.ConnectivityList
%
%   bisectN  OPTIONAL
%            Default value emtpy. A set of normals describing planes. Plane
%            k is given by all points, x, satisfying the equation
%            bisectN(k,:)*(bisectX0(k,:) - x)' = 0. The surface described 
%            by bound, is cut by the bisectors of dt AND the planes 
%            described by bisectN and bisectX0. 
%
%   bisectX0 OPTIONAL
%            Default value empty. For more information see information
%            about bisectN.
%
% RETURNS:
%   V        A mx3 array describing the vertices of the intersection of the
%            Voronoi diagram and the surface. 
%   C        A cell array. Element C{k} is an array of indices such that 
%            V(C{k}) gives the intersection of site dt.Points(k,:) and the
%            surface bound.
%   symV     A cell array describing the topology of the intersection. 
%
% EXAMPLE:
%   p = rand(5,3);
%   dt = delaunayTriangulation(p);
%   bndPts = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0,0,1; 1,0,1; 1,1,1;0,1,1];
%   bnd = delaunayTriangulation(bndPts);
%   bound.Points = bnd.Points;
%   bound.ConnectivityList = bnd.freeBoundary
%   [V,C] = clipGrid(dt,bound);
%   for i = 1:numel(C)
%     H = convhull(V(C{i},:));
%     patch('vertices',V(C{i},:),'faces',H,'facealpha',0.3)
%   end
%   axis equal
%
% SEE ALSO
%   pebi, clippedPebi3D, voronoi2mrst

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

opt = struct('bisectN', [],...
             'bisectX0',[]);
           
opt = merge_options(opt, varargin{:});

bisectN = opt.bisectN;
bisectX0 = opt.bisectX0;
assert(size(bisectN,1)==size(bisectX0,1),'# of normals must equal # of x0')

boundDim = size(bound.ConnectivityList,2);

s = dsearchn(dt.Points,sum(bound.Points(bound.ConnectivityList(1,:),:)/boundDim,1));

Q = [1, s];
V = [];
symV = cell(0);
C = cell(size(dt.Points,1),1);
CT = cell(numel(C),1);
CT{s} = 1;
E = dt.edges;
nFac  = size(bound.ConnectivityList,1);
bisId = -(nFac+(1:size(bisectN,1))');

while ~isempty(Q)
    t =  Q(end,1); s = Q(end,2);
    Q = Q(1:end-1,:);
    NC = [E(:,2)==s, E(:,1)==s];
    bisect = [find(any(NC,2)); bisId];
    
    NT = findNeighbours(bound.ConnectivityList, t);    %sum(ismember(face2vert,face2vert(t,:)),2)==2;


    n = bsxfun(@minus, dt.Points(E(NC),:), dt.Points(s,:));
    n = [bsxfun(@rdivide, n,sqrt(sum(n.^2,2))); bisectN];
    x0 =[bsxfun(@plus, dt.Points(E(NC),:), dt.Points(s,:))/2; bisectX0];

    %symT = [-find(any(bound.ConnectivityList==bound.ConnectivityList(t,1),2))'; ...
    %        -find(any(bound.ConnectivityList==bound.ConnectivityList(t,2),2))'; ...
    %        -find(any(bound.ConnectivityList==bound.ConnectivityList(t,3),2))';];
    tri2 = bound.ConnectivityList(NT,:);
    tri2 = repmat(tri2,3,1);
    tri1 = repmat(bound.ConnectivityList(t,:),9,1);
    tri1 = reshape(tri1,3,9)';
    
    IA = any(tri2==tri1,2);
    I  = mod(find(IA)-1,3)+1;
    symT = -[repmat(t,3,1), reshape(NT(I),2,3)'];
    symT = mat2cell(symT, ones(size(symT,1),1),3);
    [newVertex, symT] = clipPolygon(bound.Points(bound.ConnectivityList(t,:),:),...
                                    n,x0,symT,bisect);
    if isempty(newVertex)
      continue
    end
    symV = [symV; symT];
    C{s} = [C{s}, size(V,1)+1:size(V,1)+size(newVertex,1)];
    V = [V;newVertex];
    [Q,CT] = updateQue(Q, symT, CT, E, NC, s, t,-nFac);    
    

end

% Remove duplicate nodes.
[V, C, symV] = cleanUpGrid(V, C,symV);
end


function NT = findNeighbours(V, t)
    VT = V(t,:);
    V  = [V(1:t-1,:);nan,nan,nan;V(t+1:end,:)];
    NT = [find(sum(ismember(V,VT([1,2])),2)==2);...
          find(sum(ismember(V,VT([2,3])),2)==2);...
          find(sum(ismember(V,VT([3,1])),2)==2)];
end



function [Q, CT] = updateQue(Q, symV, CT, E, NC, s, t,tmin)
    % Find possible new cells
    symV = horzcat(symV{:})';
    
    bNew = unique(symV(symV>0));
    tNew = -unique(symV(tmin<=symV & symV<0));
    for i = 1:numel(bNew)
       if isempty(CT{E(bNew(i),NC(bNew(i),:))}) || ~any(CT{E(bNew(i),NC(bNew(i),:))}==t) %New cell facet pair
           Q = [Q; t, E(bNew(i),NC(bNew(i),:))];
           CT{E(bNew(i),NC(bNew(i),:))} = [CT{E(bNew(i),NC(bNew(i),:))}, t];
       end
    end
    for i = 1:numel(tNew)
        if ~any(CT{s}==tNew(i))
           Q = [Q; tNew(i), s];
           CT{s} = [CT{s}, tNew(i)];
        end
    end
end

function [V, C, symV] = cleanUpGrid(V, C,symV)
    % Remove duplicate vertexes
    [V,IA,IC] = uniquetol(V,1e-10,'byRows',true);
    symV = symV(IA,:);
    C = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);
    
end
