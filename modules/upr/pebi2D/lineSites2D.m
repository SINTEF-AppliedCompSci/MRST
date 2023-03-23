function [CCPts, cGs, protPts,pGs] = lineSites2D(cellConstraints, CCGridSize, varargin) 
% Places sites along given lines.
%
% SYNOPSIS:
%   [wellPts, wGs, protPts, pGs] = lineSites2D(F,faultGridSize)
%   [...] = lineSites2D(..., 'Name1', Value1, 'Name2', Value2,...)
%
% Parameters:
%   cellConstraints     - A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a well. The 
%                   values must be sorted along the line, e.g., a well 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%                   If an array has length 1 it is considered a point well.
%
%   CCGridSize    - Desired distance between sites along a well.
%                   The true distance is set such that it is a factor of 
%                   the total line length, and therefore might be slightly
%                   smaller than the desired distance.
%   
%   interpolateCC - OPTIONAL.
%                   Default value is a boolean false value, but the user
%                   can supply an individual value per line path. If false,
%                   each segment in the corresponding line curve will be
%                   represented by at least one cell. If true, the routine
%                   will interpolate along the curve, which means that
%                   sites will not necessarily fall exactly on the
%                   prescribed curve.
%
%   cfCut         - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of cell constraints. The value equals the output of the 
%                   function [~,~, cfCut] = splitAtInt2D. The value of 
%                   element i tells if the start or end site of constraint i 
%                   should be removed. If the value is 1 the end site is
%                   removed. If the value is 2 the start site is removed.
%                   If the value is 3 both the starts and end site is 
%                   removed.
%
%   wCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of cell constraints. The value equals the output of the 
%                   function [~,wCut, ~] = splitAtInt2D. The value of 
%                   element i tells if the start or end site of well i 
%                   should be removed. If the value is 1 the end site is
%                   removed. If the value is 2 the start site is removed.
%                   If the value is 3 both the starts and end site is 
%                   removed.
%
%   sePtn         - OPTIONAL.
%                   Default value is an array of zeros, Array of length
%                   equal the number of cell constraints x 2 ([s,e]). Each row gives
%                   the start and end position of the constraint interpolation.
%                   The position is relative to the steplength. If
%                   s=ones(numel(cellConstraints),1) all constraints will start their
%                   first site one step length from the start of the
%                   constrained paths. Be careful when using both wfCut and sePts
%                   as the effects add up.
%
%   protLayer     - OPTIONAL.
%                   Default set to false. If set to true a protection layer
%                   is added on both sides of the cell constraints
%                   .          .             .  Protection Layer
%                                               protD(dist between sites)
%                   .----------.-------------.  constrained path
%
%                   .          .             .  Protection Layer
% 
%  protD          - OPTIONAL.
%                   Default value CCGridSize/10. Cell array of Functions.
%                   The array should have either one function or one 
%                   function for  each cell constraint.
%                   The functions give the distance from constrained sites 
%                   and protection sites. The function is evaluated along 
%                   the constraints such that protD(0) is the start of the 
%                   constraint while protD(1) is the end of the constraint.
%
% RETURNS:
%   CCPts           Array of all generated cell constriant sites.
%   cGs             Distance between consecutive cell constraint sites.
%   protPts         Array of generated protection sites.
%   pGs             Array with the distance between CCPts(i,:) and its
%                   protection Points
%
% EXAMPLE
%   CC = {[0.2,0.2;0.8,0.8]};
%   gS = 0.1;
%   CCPts = lineSites2D(CC, gS);
%   figure(); hold on
%   plotLinePath(CC)
%   plot(CCPts(:,1), CCPts(:,2),'.r','markersize',20)
%
% SEE ALSO:
%   pebiGrid2D, compositePebiGrid2D, surfaceSites2D, pebi
%   lineSites3D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% load options
opt   = struct('cfCut',         zeros(numel(cellConstraints),1),...
               'cCut',          zeros(numel(cellConstraints),1),...
               'sePtn',         zeros(numel(cellConstraints),2),...
               'protLayer',     false, ...
               'interpolateCC', false, ...
               'protD',         {{@(p) ones(size(p,1),1)*CCGridSize/10}},...
               'CCRho',         @(x)CCGridSize*constFunc(x));
opt      = merge_options(opt,varargin{:});
cfCut    = opt.cfCut;
cCut     = opt.cCut;
sePtn    = opt.sePtn;
protD    = opt.protD;
interpCC = opt.interpolateCC;

if numel(protD) == 1
  protD = repmat(protD,numel(cellConstraints),1);num2cell(protD, 2);
end
assert(numel(protD) == numel(cellConstraints));
if numel(interpCC) == 1
  interpCC = repmat(interpCC, numel(cellConstraints),1);
end
assert(numel(interpCC)==numel(cellConstraints));

cGs     = [];
pGs     = [];
CCPts = [];
protPts = zeros(0,2);
for i = 1:numel(cellConstraints)  % create well points
  cellConstraint       = cellConstraints{i};
  
  if (size(cellConstraint,1) == 1)
      p = cellConstraint;
      wellSpace = CCGridSize;
  else
      p = interLinePath(cellConstraint, opt.CCRho, CCGridSize, sePtn(i,:), interpCC(i));
      if isempty(p)
        continue
      end
      wellSpace = sqrt(sum(diff(p,1,1).^2,2));
      wellSpace = min(wellSpace([1 1:end]), wellSpace([1:end end]));
  end
  
  keep = 1:size(p,1);
  switch cfCut(i)
    case 1
      keep = keep(1:end-1);
    case 2
      keep = keep(2:end);
    case 3
      keep = keep(2:end-1);
  end
  keepProt = keep;
  switch cCut(i)
    case 1
      keepProt = keepProt(1:end-1);
    case 2
      keepProt = keepProt(2:end);
    case 3
      keepProt = keepProt(2:end-1);
  end
  cGs     = [cGs; wellSpace(keep)];
  CCPts = [CCPts;p(keep,:)];

  if size(cellConstraint,1)>1
  if opt.protLayer
    % Calculate numerical normals
    if numel(keepProt)==1
      pK = [p(keepProt,:);cellConstraint(end,:)];
    else
      pK = p(keepProt,:);
    end
    
    newN = diff(pK);
    newN = [newN(:,2), -newN(:,1)];
    newN = bsxfun(@rdivide, newN,sqrt(sum(newN.^2,2)));
    
    if numel(keepProt)> 1
      newN = newN([1:end end],:);
    else
      pK = pK(1,:);
    end
    d = repmat(protD{i}(pK), 1,2);
    protPts = [protPts; pK + newN.*d; pK - newN.*d];
    pGs = [pGs; d(:)];
  end
  end
end

CCPts = round(CCPts*10^13)/10^13;
[CCPts,IA] = unique(CCPts,'rows'  );
cGs          = cGs(IA);

end
