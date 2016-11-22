function [wellPts, wGs, protPts,pGs] = createWellGridPoints(wellLines, wellGridSize, varargin) 
% Places well grid points along wells.
%
% SYNOPSIS:
%   [wellPts, wGs, protPts, pGs] = createWellGridPoints(F,faultGridSize)
%   [...] = createWellGridPoints(..., 'Name1', Value1, 'Name2', Value2,...)
%
% Parameters:
%   wellLines       A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a well. The 
%                   values must be sorted along the line, e.g., a well 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%                   If an array has length 1 it is considered a point well.
%
%   wellGridSize    Desired distance between well points along a well.
%                   The true distance is set such that it is a factor of 
%                   the total fault length, and therefore might be slightly
%                   smaller than the desired distance.
%   
%   wfCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of wells. The value equals the output of the 
%                   function [~,~, fwCut] = splitLines. The value of 
%                   element i tells if the start or end point of well i 
%                   should be removed. If the value is 1 the end point is
%                   removed. If the value is 2 the start point is removed.
%                   If the value is 3 both the starts and end point is 
%                   removed.
%
%   wCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of wells. The value equals the output of the 
%                   function [~,wCut, ~] = splitLines. The value of 
%                   element i tells if the start or end point of well i 
%                   should be removed. If the value is 1 the end point is
%                   removed. If the value is 2 the start point is removed.
%                   If the value is 3 both the starts and end point is 
%                   removed.
%
%   sePtn           - OPTIONAL.
%                   Default vaule array of zeros, Array of length equal the
%                   number of well x 2 ([s,e]). Each row gives the start 
%                   and end position of the well interpolation. The 
%                   position is relative the steplength. If 
%                   s=ones(numel(wellLines),1) all wells will start their 
%                   first well site one step length from the start of the
%                   well paths. Be carefull when using both wfCut and sePts
%                   as the effects adds up. 
%
%   protLayer       - OPTIONAL.
%                   Default set to false. If set to true a protection layer
%                   is added on both sides of the well
%                   .          .             .  Protection Layer
%                                               protD(dist between points)
%                   .----------.-------------.  Well path
%
%                   .          .             .  Protection Layer
% 
%  protD            - OPTIONAL.
%                   Default value wellGridSize/10. Cell array of Functions.
%                   The array should have either one function or one 
%                   function for  each well path.
%                   The functions give the distance from well sites 
%                   and protection sites. The function is evaluated along 
%                   the well path such that protD(0) is the start of the 
%                   fault while protD(1) is the end of the fault.
%
% RETURNS:
%   wellPts         Array of all generated well points.
%   wGs             Distance between consecutive well points.
%   protPts         Array of generated protection points.
%   pGs             Array with the distance between wellSite(i,:) and its
%                   protection Points
%
% EXAMPLE
%   wl = {[0.2,0.2;0.8,0.8]};
%   gS = 0.1;
%   wPts = createWellGridPoints(fl,gS);
%   figure(); hold on
%   plot(wl{1}(:,1), wl{1}(:,2))
%   plot(wPts(:,1),wPts(:,2),'.r','markersize',20)
%
% SEE ALSO:
%   pebiGrid, compositePebiGrid, createFaultGridPoints, pebi
%   createWellGridPoints3D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% load options
opt   = struct('wfCut',zeros(numel(wellLines),1),...
               'wCut',zeros(numel(wellLines),1),...
               'sePtn', zeros(numel(wellLines),2),...
               'protLayer',false, ...
               'protD', {{@(p) ones(size(p,1),1)*wellGridSize/10}},...
							 'wellRho',@(x)wellGridSize*constFunc(x));
opt   = merge_options(opt,varargin{:});
wfCut = opt.wfCut;
wCut  = opt.wCut;
sePtn = opt.sePtn;
protD = opt.protD;

if numel(protD) == 1
  protD = repmat(protD,numel(wellLines),1);num2cell(protD, 2);
end
assert(numel(protD) == numel(wellLines));

wGs     = [];
pGs     = [];
wellPts = [];
protPts = zeros(0,2);
for i = 1:numel(wellLines)  % create well points
  wellLine       = wellLines{i};
  
  if (size(wellLine,1) == 1)
      p = wellLine;
      wellSpace = wellGridSize;
  else
      p = interLinePath(wellLine, opt.wellRho,wellGridSize, sePtn(i,:));
      if isempty(p)
        continue
      end
      wellSpace = sqrt(sum(diff(p,1,1).^2,2));
      wellSpace = [wellSpace(1);wellSpace;wellSpace(end)];
      wellSpace = min([wellSpace(1:end-1), wellSpace(2:end)],[],2);
  end
  
  keep = 1:size(p,1);
  switch wfCut(i)
    case 1
      keep = keep(1:end-1);
    case 2
      keep = keep(2:end);
    case 3
      keep = keep(2:end-1);
  end
  keepProt = keep;
  switch wCut(i)
    case 1
      keepProt = keepProt(1:end-1);
    case 2
      keepProt = keepProt(2:end);
    case 3
      keepProt = keepProt(2:end-1);
  end
  wGs     = [wGs; wellSpace(keep)];
  wellPts = [wellPts;p(keep,:)];

  if size(wellLine,1)>1
  if opt.protLayer
    % Calculate numerical normals
    if numel(keepProt)==1
      pK = [p(keepProt,:);wellLine(end,:)];
    else
      pK = p(keepProt,:);
    end
    
    newN = diff(pK);
    newN = [newN(:,2), -newN(:,1)];
    newN = bsxfun(@rdivide, newN,sqrt(sum(newN.^2,2)));
    
    if numel(keepProt)> 1
      newN = [newN;newN(end,:)];
    else
      pK = pK(1,:); 
    end
    d = repmat(protD{i}(pK), 1,2);
    protPts = [protPts; pK + newN.*d; pK - newN.*d];
    pGs = [pGs; d(:)];
  end
  end
end

wellPts = round(wellPts*10^13)/10^13;
[wellPts,IA] = unique(wellPts,'rows'  );
wGs          = wGs(IA);

end
