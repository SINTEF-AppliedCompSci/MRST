function [p, symP] = clipPolygon(p, n, x0, symP, bisector,varargin)
% Clip a polygon against a set of bounding planes

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
TOL = 100*eps;

assert(size(n,1)==size(x0,1),'unconsistent size of n and x0');
assert(size(p,1)>2, 'A polygon needs more than 2 vertexes');

opt = struct('noSym',false);
opt = merge_options(opt,varargin{:});

%IC = zeros(size(p,1),1);
% symP = [-4,-2,-1;
%         -3,-2,-1;
%         -4,-3,-1];
for i = 1:size(n,1)
    
    % find distance from point to plane
    d = bsxfun(@ minus, p, x0(i,:))*n(i,:)';
    
    if all(d>TOL)
        p  = [];
        %IC = [];
        symP = [];
        return
    elseif all(d<TOL)
      continue
    end
    p = [p(end,:);p];
    d = [d(end);  d];
    symP = [symP(end); symP];

    c1 = d(1:end-1)<TOL & d(2:end)<TOL;
    c2 = sign(d(1:end-1))~=sign(d(2:end));
    c3 = (c1) | (c2 & d(2:end)<-TOL);
    c4 = false(size(p,1)-1,1);
    
    cVert = abs(d)<TOL;
    atStr = cVert(1:end-1);
    atEnd = cVert(2:end);

    c3(atEnd) = false;
    c4(atEnd) = true;
    c2(atEnd) = false;
    c2(atStr) = false;
    
    c2 = find(c2);
    c3 = find(c3);    
    c4 = find(c4);
    shift = d(c4)>TOL;
    c4 = mod(c4+shift-1,size(p,1)-1)+1;
    if ~opt.noSym
      symP = [intersectUnion(symP,c2, bisector(i));         ...
              intersectUnion(symP,c4, bisector(i));        ...
              symP(c3+1,:)];
    end

    alpha = abs(d(c2))./(abs(d(c2))+abs(d(c2+1)));
    if isempty(alpha)
      p = [p(c4+~shift,:);p(c3+1,:)]; 
    else
      p  = [bsxfun(@times,p(c2+1,:)-p(c2,:),alpha)+p(c2,:); ...
            p(c4+~shift,:);                                 ...
            p(c3+1,:)];
            
    end
    %IC = [i*ones(size(c2,1),1);IC(c3)];
    c = [c2;c4;c3];
    [~,I] = sort(c);
    p   = p(I,:);
    %IC = IC(I);
    if ~opt.noSym
      symP = symP(I,:);
    end
end


 
end


function [B] = intersectUnion(A, idx, b)
    B = cell(size(idx,1),1);
    for i = 1:size(idx,1)
       B{i,:} = [intersect(A{idx(i)},A{idx(i)+1}),b]; 
    end
    
end
