function [V, C, symV] = clipGrid(dt,bound,varargin)

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


function [symV] = updateSym(localSym, NC, NT)
    if localSym<0
        symV = -NT(-localSym);
    else
        symV = NC(localSym);
    end
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
%     symV = cellfun(@sort, symV,'UniformOutput', false);
%     symV = cellfun(@(c) num2str(c'), symV,'UniformOutput', false);
%     [~,IA,IC] = unique(symV);
%     VO = V;
%     CO = C;
%     C = cellfun(@(c) unique(IC(c)'), C,'UniformOutput',false);
%     V = V(IA,:);
    [V,IA,IC] = uniquetol(V,1e-10,'byRows',true);
    symV = symV(IA,:);
    C = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);
    
end
