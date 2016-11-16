function [F] = createFaultGridPoints(faultLines,faultGridSize, varargin) 
% Places fault grid points on both sides of faults.
%
% SYNOPSIS:
%   F = createFaultGridPoints(F,faultGridSize)
%   F = createFaultGridPoints(..., 'Name1', Value1, 'Name2', Value2, ...)
%
% Parameters:
%   faultLine       A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a fault. The 
%                   values must be sorted along the line, e.g., a fault 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%   faultGridSize   Desired distance between fault points along a fault.
%                   The true distance is set such that it is a factor of 
%                   the total fault length, and therefore might be slightly
%                   smaller than the desired distance.
%   
%   circleFactor    - OPTIONAL.
%                   Default value 0.6. The ratio between faultGridSize and 
%                   the radius of the circles used to create the points. If
%                   circleFactor is increased, the distance between the 
%                   fault and the points are increased. Valid values for
%                   circleFactor are in the interval (0.5,1.0).
%
%   distFun         - OPTIONAL.
%                   Default value @(x) faultGridSize*constFunc(x). Function
%                   handle that specify a relative fault grid size. If
%                   distFun = @(x) 1-0.5*x(:,1), in the unit square, then
%                   any any faults on the right side will have about twice
%                   as many grid points as equivalent faults placed on the
%                   left side. 
%
%   fCut            - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of faults. The value equals the output of the 
%                   function [~, fCut,~]=splitFaults. The value of element
%                   i tells if fault i share a start or end point with 
%                   other faults. If the value is 1 it share an end point.
%                   If the value is 2 it share a start point. If the value
%                   is 3 it share both a start point and an end point.
%
%   fwCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of faults. The value equals the output of the 
%                   function [~,~, fwCut] = splitFaults. The value of 
%                   element i tells if the start or end point of fault i 
%                   should start half a step length from the start/end. If
%                   the value is 1 it starts half a step length from the 
%                   end. If the value is 2 it starts half a step length 
%                   from the start. If the value is 3 it starts half a step
%                   length from both the start and end.
%
%   mergeTol        - OPTIONAL.
%                   Default value zero. Relative allowed distance between
%                   two circles on different fault paths. If the relative
%                   distance is less than a 
%                   (c1 - c2).^2 <(mergeTol*fh(c1))^2 the two circles c1 
%                   and c2 are merged.
% RETURNS:
%   F               Struct with elements:
%     F.f.pts       Point coordinates
%     F.f.Gs        Grid spacing for each fault point. This is the distance
%                   between the two points on oposite sides of a fault.
%     F.f.c         Map from fault points to circles
%     F.f.cPos      Fault point i was created using circle
%                   F.f.c(F.f.cPos(i):F.f.cPos(i+1)-1).
%     F.c.CC        Coordinates to circle centers
%     F.c.R         Radius of circles
%     F.c.f         Map from circles to fault points
%     F.c.fPos      Circle i created fault points
%                   F.c.f(F.c.fPos(i):F.c.fPos(i+1)-1)
%     F.c.l         Map from circles to fault
%     F.c.lPos      Circle i lies on fault paths
%                   F.l.c(F.l.cPos(i):F.l.cPos(i+1)-1) 
%     F.l.f         Map from fault lines to fault points
%     F.l.fPos      Fault points F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1) is 
%                   created for fault i.
%     F.l.c         Map from fault paths to circles
%     F.l.cPos      Circle F.l.c(F.l.cPos(i):F.l.cPos(i+1)-1) is created
%                   for fault i.
%     F.l.l         faultLine from input.
%     
%
% EXAMPLE
%   fl = {[0.2,0.2;0.8,0.8]};
%   gS = 0.1;
%   F = createFaultGridPoints(fl,gS);
%   figure(); hold on
%   plot(fl{1}(:,1), fl{1}(:,2))
%   plot(F.f.pts(:,1),F.f.pts(:,2),'.r','markersize',20)
%
% SEE ALSO:
%   pebiGrid, compositePebiGrid, createWellGridPoints, splitFaults, pebi.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Load options

fh  = @(x) faultGridSize*constFunc(x);
opt = struct('distFun',     fh,  ...
             'circleFactor', 0.6, ...
             'fCut',         zeros(numel(faultLines),1), ...
             'fwCut',        zeros(numel(faultLines),1), ...
             'mergeTol',     0);
opt = merge_options(opt, varargin{:});

fh           = opt.distFun;
circleFactor = opt.circleFactor;
fCut         = opt.fCut;
fwCut        = opt.fwCut;


% Initialize variables.
F.f.Gs    = [];                % Fault point grid size
F.f.pts   = [];                % Fault points
F.f.c     = [];                % Map from a fault to the circle center
F.f.cPos  = 1;

F.c.CC    = [];                % Center of circle used to create fault pts
F.c.R     = [];                % Radius of the circle
F.c.f     = [];                % Map from the circle center to a fault  
F.c.fPos  = 1;
F.c.l     = [];
F.c.lPos  = 1;

F.l.f     = [];
F.l.fPos  = 1;                  % Map frm fault lines to fault points
F.l.l     = faultLines;
F.l.nFault = numel(faultLines);



for i = 1:F.l.nFault
  faultLine     = F.l.l{i};
  sePtn         = .5*[fwCut(i)==2|fwCut(i)==3; fwCut(i)==1|fwCut(i)==3];
  [p, fracSpace, fCi, fRi, f2ci,cPos, c2fi,fPos] =   ...
                          faultLinePts(faultLine,     ... 
                                       faultGridSize,...
                                       circleFactor, ...
                                       fCut(i),sePtn,fh);
  nl = size(p,1)/2;
  if nl==0 % No fault points created
    F.l.fPos = [F.l.fPos; F.l.fPos(end)];
    continue
  end
  F.f.Gs   = [F.f.Gs;fracSpace];
  F.l.fPos = [F.l.fPos; size(F.f.pts,1)+1+size(p,1)];
  F.f.pts  = [F.f.pts;p];
  F.c.CC   = [F.c.CC; fCi];
  F.c.R    = [F.c.R; fRi]; 
  F.f.c    = [F.f.c; f2ci + size(F.c.fPos,1)-1];
  F.f.cPos = [F.f.cPos; cPos(2:end) + F.f.cPos(end)-1];
  F.c.f    = [F.c.f; c2fi + size(F.f.pts,1)-nl*2];
  F.c.fPos = [F.c.fPos; fPos(2:end) + F.c.fPos(end)-1];
  F.c.l    = [F.c.l;repmat(i,size(fCi,1),1)];
end
F.l.f = (1:F.l.fPos(end)-1)';
F.c.lPos = (1:size(F.c.CC,1)+1)';

% Add well-fault intersections
if ~isempty(F.c.CC)
  cutAtEnd = fwCut==1|fwCut==3;
  cutAtStr = fwCut==2|fwCut==3;
    
  endCirc = F.f.c(F.f.cPos(F.l.fPos([false;cutAtEnd]))-1); %F.l.f is has not changed
  strCirc = F.f.c(F.f.cPos(F.l.fPos(cutAtStr)));

  p       = circCircInt(F.c.CC(strCirc,:), F.c.R(strCirc),...
                        F.c.CC(endCirc,:), F.c.R(endCirc));
  l2fId1  = repmat(F.l.fPos(cutAtStr),1,2);
  l2fId2  = repmat(F.l.fPos([false;cutAtEnd]),1,2); %Not -1 because of how insertVec works
  l2fId   = reshape([l2fId1,l2fId2]',[],1);
  fId     = (size(F.f.pts,1)+1:size(F.f.pts,1) + size(p,1))';
  fId     = reshape(fId,2,[]);
  fId     = repmat(fId,2,1);
  cId     = reshape([strCirc,endCirc]',[],1);
  c2fId   = repmat(F.c.fPos(cId)',2,1);
  c2fId   = c2fId(:);

  F.l.f   = insertVec(F.l.f, fId(:), l2fId);
  F.f.pts = [F.f.pts;p];
  nGs     = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
  F.f.Gs  = [F.f.Gs;reshape(nGs(:,1:2:end),[],1)];
  F.c.fPos= F.c.fPos + ...
    cumsum(accumarray([cId+1;size(F.c.fPos,1)],2*[ones(1,size(cId,1)),0]));
  F.c.f   = insertVec(F.c.f, fId(:), c2fId);
  cId     = repmat(reshape(cId,2,[]),2,1);
  F.f.cPos= [F.f.cPos; F.f.cPos(end)+2*cumsum(ones(size(p,1),1))];
  F.f.c   = [F.f.c;cId(:)];
  
  F.l.fPos(2:end) = F.l.fPos(2:end) + 2*cumsum(cutAtEnd+cutAtStr);
end


% Merge fault intersections
if ~isempty(F.f.pts)
  % Remove duplicate fault centers
  F.c.CC = floor(F.c.CC*10^15)/10^15;
  [~, IA, IC] = unique(F.c.CC,'rows');
  F.c.CC      = F.c.CC(IA,:);
  F.c.R       = F.c.R(IA);
  [~,I]       = sort(IC);
  from2f      = [F.c.fPos(1:end-1), F.c.fPos(2:end)-1];
  from2f      = from2f(I,:);
  F.c.f       = F.c.f(mcolon(from2f(:,1),from2f(:,2)));
  from2l      = [F.c.lPos(1:end-1), F.c.lPos(2:end)-1];
  from2l      = from2l(I,:);
  F.c.l       = F.c.l(mcolon(from2l(:,1),from2l(:,2)));
  fNum        = diff(F.c.fPos);
  lNum        = diff(F.c.lPos);
  F.c.fPos    = cumsum([1; accumarray(IC,fNum)]);
  F.c.lPos    = cumsum([1; accumarray(IC,lNum)]);
  F.f.c       = IC(F.f.c);
  
  % Merge intersections
  [F]         = fixIntersections(F, fh, circleFactor,opt.mergeTol);
end
end

function [Pts, gridSpacing, circCenter, circRadius, f2c,f2cPos, c2f,c2fPos] = ...
    faultLinePts(faultLine, fracDs, circleFactor, isCut,sePtn, fh) 

    assert(0.5<circleFactor && circleFactor<1)
    assert(size(faultLine,1)>1 && size(faultLine,2)==2);
    
    % interpolate fault line to get circle centers. 
    circCenter = interLinePath(faultLine, fh, fracDs,sePtn);
    numOfFracPts = size(circCenter,1)-1;
    
    % Test if faultLine is to short
    if numOfFracPts == 1
      d = sqrt(sum((circCenter(2,:)-circCenter(1,:)).^2, 2));
      if d < 0.8*fh((circCenter(2,:)+circCenter(1,:))/2);
        Pts         = [];
        gridSpacing = [];
        circCenter  = [];
        circRadius  = [];
        f2c         = [];
        f2cPos      = [];
        c2f         = [];
        c2fPos      = [];
        return
      end
    end
    % Calculate the line length and circle radii. If you experience
    % imaginary faultOffset you might want to try the max lineLength
    % instead of the mean.
    lineLength = sqrt(sum((circCenter(2:end,:)-circCenter(1:end-1,:)).^2, 2));
    circRadius = circleFactor*[lineLength(1); ...
                              (lineLength(1:end-1) + lineLength(2:end))/2; ...
                               lineLength(end)];
                             
    switch isCut
      case 1
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
      case 2
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
      case 3
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
    end
    
    % Calculate the crossing of the circles
    bisectPnt = (lineLength.^2 - circRadius(2:end).^2 + circRadius(1:end-1).^2)...
                ./(2*lineLength);
    faultOffset = sqrt(circRadius(1:end-1).^2 - bisectPnt.^2);
    n1 = (circCenter(2:end,:)-circCenter(1:end-1,:))./repmat(lineLength,1,2); %Unit vector
    n2 = [-n1(:, 2), n1(:,1)];                                                %Unit normal
    
    % Set fault points on left and right side of fault
    left   = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             + bsxfun(@times, faultOffset, n2);
    right  = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             - bsxfun(@times, faultOffset, n2);
         
    % Put together result
    Pts = [right;left];
    f2c = cumsum(accumarray((1:2:size(Pts,1)+1)',1));
    f2c = repmat(f2c(2:end),2,1); 
    f2cPos = (1:2:numel(f2c)+1)';
    nf  = size(left,1);
    c2f = [   nan,          nan,         1,       nf+1;...
           (1:nf-1)', (nf+1:2*nf-1)', (2:nf)', (nf+2:2*nf)';...
              nf,           2*nf,       nan,     nan]';
    c2f = c2f(3:end-2)';
    c2fPos = [1;(3:4:numel(c2f))';numel(c2f)+1];
    gridSpacing = 2*[faultOffset;faultOffset];
end



function [F] = fixIntersections(F,fh, circFac, mergeTol)
  if ~all(diff(F.f.cPos)==2)
    warning('Can only merge sites that are created from exactly 2 circles');
    return
  end
  % Find conflict circles
  I = conflictCircles(F.f.pts, F.c.CC, F.c.R);
  
  circ    = find(~cellfun(@isempty,I));
  if isempty(circ)
    return
  end

  circNum = cellfun(@numel, I(circ));
  id      = zeros(sum(circNum),1);
  circPos = cumsum([1; circNum]);
  id(circPos(1:end-1)) = 1;
  id      = cumsum(id);
  circ    = [circ(id), F.f.c([F.f.cPos(vertcat(I{circ})),...
             F.f.cPos(vertcat(I{circ}))+1])];
  
  % Find shared circle
  [neigh,neighPos] = findNeighbors(circ(:,1),F);%F.c.f,F.c.fPos, F.f.c,F.f.cPos);
  if any(~all(diff(neighPos)==2))
    warning('Something went wrong. Can not guarantee conformity')
    remC = diff(neighPos)~=2;
    circ(remC,:) = [];
    remN = mcolon(neighPos([remC;false]), neighPos([false;remC])-1);
    neigh(remN) = [];
    neighPos = neighPos - cumsum([0;remC]);
    neighPos(remC) = [];
    assert(all(diff(neighPos)==2) )
  end
  
  neigh  = reshape(neigh,2,[])';
  shared = 2*any(bsxfun(@eq, neigh, circ(:,2)),2) ...
          +3*any(bsxfun(@eq, neigh, circ(:,3)),2);
  keep   = find(shared);
  circ   = circ(keep,:);
  swap   = shared(keep)==2;                % Set shared circle at third row
  circ(swap,:) = [circ(swap,1),circ(swap,3),circ(swap,2)];
  % circ(:,1:2) is now conflict circles pairs which we either shrink the 
  % radii to or merge. circ(:,3) is the conflic circles common neighbor.
  
  % Remove duplicate pairs
  [~,IA] = unique(sort(circ,2),'rows');
  circ   = circ(IA,:);
  
  % Calculate new radii
  line = [F.c.CC(circ(:,3),:),reshape(mean(reshape(F.c.CC(circ(:,1:2)',:),2,[]),1),[],2)];
  int  = lineCircInt(F.c.CC(circ(:,3),:),F.c.R(circ(:,3)), line);

  merge = (F.c.CC(circ(:,1)) - F.c.CC(circ(:,2))).^2 ...
            <(mergeTol*fh(F.c.CC(circ(:,1),:))).^2;
  if any(merge)
    
    mC = [circ(:,1), circ(:,2)];
    F = mergeCirc(F,mC,fh,circFac);

    F = fixIntersections(F,fh, circFac);

    return
  end
  % set radius to smallest
  R = sqrt(sum((F.c.CC(circ(:,1:2),:)-[int;int]).^2,2));
  Rold = F.c.R;
  for i = 1:numel(R)
    if R(i)<F.c.R(circ(i))
      F.c.R(circ(i)) = R(i);
    end
  end
  c = unique(circ(:,1:2));
  if size(c,2)>1
    c = c';
  end
  % Calculate new Pts
  map = arrayfun(@colon,F.c.fPos(c),F.c.fPos(c+1)-1,'uniformOutput',false)';
  map = horzcat(map{:})';
  fId = F.c.f(map);
%   fId = unique(fId(:));
%   fId = fId(~isnan(fId));

  [neigh,neighPos] = findNeighbors(c, F);%F.c.f,F.c.fPos, F.f.c,F.f.cPos);
  assert(all(diff(neighPos)==2));
  neigh = reshape(neigh,2,[])';
  
  p = circCircInt(F.c.CC(c,:), F.c.R(c),...
                 reshape(F.c.CC(neigh',:)',4,[])',reshape(F.c.R(neigh),[],2));
  reN = any(p - conj(p)==0,2);
  p = p(reN,:);
  fId = fId(reN);

  if any(~reN)
    F.c.R = Rold;
    circNum = rldecode(1:size(F.c.CC,1), diff(F.c.fPos),2);
    mC = circNum(map(~reN));  
    %mC = unique(mC, 'stable');
    [mC,~,IC] = intersect(mC, circ(:,1:2));
    m = size(circ,1);
    mC = [mC, circ(mod(IC+m-1,2*m)+1)]; %c(IC + (-1+2*bitget(IC,1)))];    
    if size(mC,1)==1
      mC = reshape(mC,[],2);
    end
    mC = unique(sort(mC,2),'rows');

    F = mergeCirc(F,mC,fh,circFac);
    F = fixIntersections(F,fh, circFac,mergeTol);

    return
  end
  p = real(p);
  F.f.pts(fId,:) = p;  
  nGs            = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
  F.f.Gs(fId)    = reshape(nGs(:,1:2:end),[],1);
  map            = [F.f.cPos(fId),F.f.cPos(fId)+1]';
  
  F.f.c(map(:))  = [c';neigh(:,1)';c';neigh(:,1)';c';neigh(:,2)';c';neigh(:,2)'];%reshape(repmat([c',c';neigh(:,1)',neigh(:,2)'],2,1),2,[]);
  p = round(F.f.pts*10^13)/10^13;
  [~, IA, IC] = unique(p,'rows');
  F.f.pts = F.f.pts(IA,:);
  F.f.Gs = F.f.Gs(IA);
  [~,I] = sort(IC);
  
  mapc = [F.f.cPos(1:end-1), F.f.cPos(2:end)-1];
  mapc = mapc(I,:);
  mapc = arrayfun(@colon, mapc(:,1),mapc(:,2),'uniformOutput',false);
  F.f.c = F.f.c(cell2mat(mapc'));
  cNum = diff(F.f.cPos);
  F.f.cPos = cumsum([1;accumarray(IC,cNum)]);
  F.c.f = IC(F.c.f);
  
  F.l.f = IC(F.l.f);
  end


function [F] = mergeCirc(F,c,fh,circFac)
  rc = reshape(c', [],1);
  C1 = c(:,1);
  C2 = c(:,2);
  newCC = (F.c.CC(C1,:) + F.c.CC(C2,:))/2;
  newR = 1.1*circFac*fh(F.c.CC(C1,:));
  [neigh,neighPos] = findNeighbors(rc,F);
  assert(all(diff(neighPos)==2));

  neigh  = reshape(neigh,4,[])';
  N1 = neigh(:,1:2);
  N2 = neigh(:,3:4);
  
  n  = bsxfun(@ne, N1(:,1), N2) & bsxfun(@ne, N1(:,2), N2);
  N2 = N2';  
  N  = [N1(:,1:2), N2(n')];

  newPts = circCircInt(newCC, newR, reshape(F.c.CC(N',:)',6,[])',reshape(F.c.R(N),[],3));
  
  c1Pos = [F.c.lPos(C1),F.c.lPos(C1+1)-1];
  c2Pos = [F.c.lPos(C2),F.c.lPos(C2+1)-1];
  c1l   = F.c.l(mcolon(c1Pos(:,1), c1Pos(:,2)));
  c2l   = F.c.l(mcolon(c2Pos(:,1), c2Pos(:,2)));
  c2S   = diff(c2Pos,1,2)+1;
  c1S   = diff(c2Pos,1,2)+1;
  clNew = insertVec(c1l, c2l,repelem((1:size(c2Pos,1))',c2S));
  lPosN = cumsum([1;c2S+c1S]);
  
  % Remove circles
  F.c.CC(rc,:) = [];
  F.c.R(rc) = [];

  mapf = mcolon(F.c.fPos(rc),F.c.fPos(rc+1)-1);
  mapl = mcolon(F.c.lPos(rc),F.c.lPos(rc+1)-1);
  rf = F.c.f(mapf);
  LIA = ismember(F.c.f,rf);
  F.c.f(LIA) = [];
  F.c.l(mapl) = [];
  
  lSize = diff([F.c.lPos(rc), F.c.lPos(rc+1)],1,2);  
  I1 = zeros(size(F.c.lPos,1)-1,1);
  I1(rc) = -lSize;
  I1 = [0;cumsum(I1)];

  I2 = zeros(size(F.c.fPos,1)-1,1);
  I2(rc) = -1;
  I2     = [0; cumsum(I2)];
 
  if size(N,1)==1
    N = N + I2(N)';
  else
    N = N + I2(N);
  end
  F.f.c = F.f.c + I2(F.f.c);
  % Remove sites
  F.f.pts(rf,:) = [];
  F.f.Gs(rf,:) = [];
  mapc = mcolon(F.f.cPos(rf),F.f.cPos(rf+1)-1);
  F.f.c(mapc) = [];
  
  I4 = zeros(size(F.f.cPos,1)-1,1);
  I4(rf) = -1;
  I4 = [0;cumsum(I4)];
  
  F.c.f = F.c.f + I4(F.c.f);
  
  [~,I] = sort(F.l.f);
  id = (1:numel(F.l.f))';
  IC = zeros(numel(F.l.f),1);
  IC(I) = id;
  LIA = ismember(F.l.f(I),rf);
  sub = cumsum(LIA);
  sub = sub(IC);
  LIA = LIA(IC);
  
  F.l.f = F.l.f - sub;
  
  F.l.f(LIA) = [];
  sub2 = cumsum(LIA);
  F.l.fPos(2:end) = F.l.fPos(2:end) - sub2(F.l.fPos(2:end)-1);
  
  F.c.lPos = F.c.lPos + I1;
  F.c.lPos(rc+1) = [];

  F.f.cPos = cumsum([1;accumarray(F.c.f,1)]);
  F.c.fPos = cumsum([1;accumarray(F.f.c,1)]);
  
  % Add circles
  F.c.CC = [F.c.CC; newCC];
  F.c.R  = [F.c.R; newR];
  % Add sites
  F.f.pts = [F.f.pts; newPts];
  nGs     = repmat(sqrt(sum(diff(newPts).^2,2)),1,2)';
  F.f.Gs  = [F.f.Gs;reshape(nGs(:,1:2:end),[],1)];
  % Add mappings
  F.f.cPos = [F.f.cPos; F.f.cPos(end) + (2:2:2*size(newPts,1))'];
  c = (size(F.c.CC,1)-size(newCC,1)+1:size(F.c.CC,1))';
  f2c = reshape([c';N(:,1)';c';N(:,1)'; ...
                 c';N(:,2)';c';N(:,2)'; ...
                 c';N(:,3)';c';N(:,3)'],2,[]);
  F.f.c = [F.f.c; f2c(:)];
  
  f = (size(F.f.pts,1)-size(newPts,1)+1: size(F.f.pts,1))';
  F.c.f = [F.c.f; f];

  cf = reshape(N',[],1);
  cf = repmat(cf,1,2);
  cf = reshape(cf',[],1);
  [cf,I] = sort(cf); % we need to sort; if circ 6 is empty fPos(cf) both fPos(6) and fPos(8) maps to the same position
  f = f(I);
  F.c.f = insertVec(F.c.f,f,F.c.fPos(cf));
  
  F.f.cPos = cumsum([1;accumarray(F.c.f,1)]);
  F.c.fPos = cumsum([1;accumarray(F.f.c,1)]);
  
  F.c.l = [F.c.l;clNew];
  F.c.lPos = [F.c.lPos; F.c.lPos(end)+lPosN(2:end)-1]; 
end


function [neigh,neighPos] = findNeighbors(c, F)
% Find the neighbor circles of c.
c2f = F.c.f;
c2fPos = F.c.fPos;
f2c = F.f.c;
f2cPos = F.f.cPos;
map   = arrayfun(@colon, c2fPos(c),c2fPos(c+1)-1,'uniformoutput',false);
pId   = cellfun(@(c) c2f(c), map,'uniformOutput',false);
neighMap = cellfun(@(c) cell2mat(arrayfun(@colon, f2cPos(c),f2cPos(c+1)-1,'uniformOutput',false)')...
                    ,pId,'uniformOutput',false);
neigh = cellfun(@(c) f2c(c),neighMap,'uniformOutput',false);
neigh = cellfun(@unique, neigh,'uniformOutput',false);
neigh = arrayfun(@(i) neigh{i}(neigh{i}~=c(i)),1:numel(neigh),'uniformOutput',false)';
neighPos = cumsum([1;cellfun(@numel, neigh)]);
neigh = vertcat(neigh{:});
end

function [p] = circCircInt(CC1, CR1, CC2,CR2)
if isempty(CC1) || isempty(CC2)
  p = [];
  return
end
% Expand matrices for computation
CC1 = repmat(CC1, 1,size(CC2,2)/size(CC1,2));
CR1 = repmat(CR1, 1,size(CR2,2)/size(CR1,2));
CC1 = reshape(CC1',2,[])';
CC2 = reshape(CC2',2,[])';
CR1 = reshape(CR1',1,[])';
CR2 = reshape(CR2',1,[])';

d = sqrt(sum((CC1 - CC2).^2,2));              % Distance between centers
bisectPnt = (d.^2 - CR2.^2 + CR1.^2)./(2*d);  % Mid-Point
faultOffset = sqrt(CR1.^2 - bisectPnt.^2);    % Pythagoras
n1 = (CC2-CC1)./repmat(d,1,2);                % Unit vector
n2 = [-n1(:, 2), n1(:,1)];                    % Unit normal

% Set right left and right intersection points
left   = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         + bsxfun(@times, faultOffset, n2);
right  = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         - bsxfun(@times, faultOffset, n2);

% Put result together
p = reshape([right,left]',2,[])';

end

function [p] = lineCircInt(CC, CR, line)
vec = line(:,3:4) - line(:,1:2);
c2l = line(:,1:2) - CC;
a   = dot(vec,vec,2);
b   = 2*dot(c2l,vec,2);
c   = dot(c2l,c2l,2) - dot(CR,CR,2);

dist    = (b.*b - 4*a.*c);
lineHit = dist>=0;
distSqr = sqrt(dist(lineHit));

%t(:,1) = -b - distSqr./(2*a); % This is the intersection on wrong side.
t = -b + distSqr./(2*a);
p = bsxfun(@times,vec,t) + line(:,1:2);
end


