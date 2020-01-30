function [splitL1, L1Cut, L1L2Cut, IC] = splitAtInt2D(L1, L2)
% Split paths at all intersections. 
%
% SYNOPSIS:
%   [splitL1, fCut, fwCut] = splitAtInt2D(L1, L2)
%
% PARAMETERS;
%   L1              A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear path. Each path in L1 is split at
%                   all intersections with other paths in L1 and all paths
%                   in L2. The values must be sorted along the path, e.g., 
%                   a  path consisting of two lines 
%                   would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%   L2              A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a path.  The 
%                   values must be sorted along the line, e.g., a cell constraint 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%                   If an array has length 1 it is considered a point.
% RETURNS:
%   splitL1         A cell of arrays. The arrays are the cut paths and 
%                   does not contain any intersections with L1 or L2,
%                   except possible at the start and/or end points.
%
%   L1Cut           Array of length equal splitL1. The value of element
%                   i tells if the returned path i has an intersection 
%                   with L1. If the value is 0, the fault has no 
%                   intersections. If the value is 1 it has an intersection 
%                   at the end point. If the value is 2 it has an 
%                   intersection at the start point. If the value is 3 it 
%                   has an intersection at both the start and end point.
%
%   L1L2Cut         Array of length equal splitL1. The value of element
%                   i tells if the returned path i has an intersection 
%                   with L2. If the value is 0, the fault has no 
%                   intersections. If the value is 1 it has an intersection 
%                   at the end point. If the value is 2 it has an 
%                   intersection at the start point. If the value is 3 it 
%                   has an intersection at both the start and end point.
%
% EXAMPLE
%   L1 = {[0.2,0.2;0.8,0.8], [0.2,0.5;0.8,0.5]};
%   L2 = {[0.3,0.8;0.3,0.2]};
%   gS = 0.1;
%   [L1,L1Cut,L1L2Cut] = splitAtInt2D(L1,L2);
%   figure(); hold on
%   name = {};
%   for i = 1:numel(L1)
%     plot(L1{i}(:,1), L1{i}(:,2))
%   end
%   plotLinePath(L2)
%   disp(L1Cut)
%   disp(L1L2Cut)
%
% SEE ALSO:
%   splitWells, surfaceSites2D, compositePebiGrid2D, pebiGrid2D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
if isempty(L1)
  splitL1 = L1;
  [L1Cut, L1L2Cut, IC] = deal([]);
  return
end
splitL1 = cell(0);
tmpCut  = cell(0);
IC      = (1:numel(L1))';
for i = 1:numel(L1)
  [splitL1{i},~,tmpCut{i}] = splitLines(L1(i),L1([1:i-1,i+1:end]));  
end
n       = cellfun(@numel,splitL1);
IC      = repelem(IC,n);
splitL1 = horzcat(splitL1{:});
tmpCut  = vertcat(tmpCut{:});
[splitL1,cutId,L1L2Cut,IC2] = splitLines(splitL1, L2);
IC  = IC(IC2);

L1Cut = zeros(sum(sum(cutId, 1) + 2),1); % initialize array
numEl = 0;
for i = 1:size(tmpCut,1)
  switch tmpCut(i)
    case 0
      numNew = sum(cutId(i,:)) + 1;
      toNum  = numEl + numNew;
      L1Cut(numEl+1:toNum,1) = [zeros(numNew-1,1); 0];
    case 1
      numNew = sum(cutId(i,:)) + 1;
      toNum  = numEl + numNew;
      L1Cut(numEl+1:toNum,1) = [zeros(numNew-1,1); 1];
    case 2
      numNew = sum(cutId(i,:)) + 1;
      toNum  = numEl + numNew;
      L1Cut(numEl+1:toNum,1) = [2; zeros(numNew-1,1)];
    case 3 
      if sum(cutId(i,:))
        numNew = sum(cutId(i,:)) + 1;
        toNum  = numEl + numNew;
        L1Cut(numEl+1:toNum,1) = [2; zeros(numNew-2,1); 1];
      else
        toNum = numEl + 1;
        L1Cut(numEl+1,1) = 3;
      end
  end
  numEl = toNum;
end
L1Cut = L1Cut(1:numEl); % Remove values not assigned.


%% Remove lines that are one point
for i = 1:numel(splitL1)
  [~, IA,~] = unique(round(splitL1{i}*10^13)/10^13,'rows','stable');
  splitL1{i} = splitL1{i}(IA, :);
end
if numel(splitL1)>0
  numPts  = cellfun(@(c) size(c,1),splitL1);  
  keep    = numPts>1;
  splitL1 = splitL1(keep);
  L1Cut   = L1Cut(keep);
  L1L2Cut = L1L2Cut(keep);
  IC      = IC(keep);
end
