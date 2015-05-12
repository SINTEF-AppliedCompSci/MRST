function CG = addCoarseCenterPointsFrac(CG,varargin)
% addCoarseCenterPointsFrac(CG) adds coarse nodes to the fracture coarse
% grid CG
%
% SYNOPSIS:
%   CG = addCoarseCenterPointsFrac(CG)
%   CG = addCoarseCenterPointsFrac(CG, 'pn1', 'pv1', ...)
%
% REQUIRED PARAMETERS:
%
%   CG  - Fracture coarse grid (supplied by 'generateCoarseGrid') with
%         geometry information (computed through 'coarsenGeometry').
%
% OPTIONAL PARAMETERS (supplied in 'key'/'value' pairs ('pn'/'pv' ...)):
%
%   option - option specifies the type of coarse grid points to use for
%            computing the coarse node location. Possible self-explanatory
%            values are: 
%
%            (a) 'useCoarseFaceCentroids'
%            (b) 'useCoarseCellCentroids'
%            (c) 'useFineCellCentroids'
%               
%            passed as character arrays.
%
%   meantype - type of mean, of the points specified by option, to use for
%              computing the coarse node location. Valid values include
%              'geometric' and 'arithmetic' passed as character arrays.
%
% RETURNS:
%   CG - Fracture coarse grid with CG.cells.centers added as a new field.
%
% SEE ALSO:
%   storeInteractionRegionFrac

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


Gf = CG.parent;
opt = struct('option', 'useCoarseFaceCentroids','meantype','geometric'); 
opt = merge_options(opt, varargin{:});
if strcmp(opt.option, 'useCoarseFaceCentroids')
    coarse_faces = rldecode(1 : CG.cells.num, diff(CG.cells.facePos), 2) .';
end
centers = zeros(CG.cells.num,1);
for i = 1:CG.cells.num
    search = find(CG.partition==i);
    if strcmp(opt.option, 'useCoarseFaceCentroids')
        faces = CG.cells.faces(coarse_faces(:,1)==i,1);
        pts = CG.faces.centroids(faces,:);
        if strcmp(opt.meantype,'geometric')
            gm = geometricMedian(pts);
        else
            gm = mean(pts);
        end
        [cx,cy] = deal(gm(1),gm(2));
    elseif strcmp(opt.option, 'useCoarseCellCentroids')
        [cx,cy] = deal(CG.cells.centroids(i,1),CG.cells.centroids(i,2));
    elseif strcmp(opt.option, 'useFineCellCentroids')
        x = Gf.cells.centroids(search,1);
        y = Gf.cells.centroids(search,2); pts = [x,y];
        if strcmp(opt.meantype,'geometric')
            gm = geometricMedian(pts);
        else
            gm = mean(pts);
        end
        [cx,cy] = deal(gm(1),gm(2));
    else
        error(['Wrong value for key ''option'' in argument list. Can be either one of ',...
            '''useCoarseFaceCentroids'', ''useCoarseCellCentroids'' or ''useFineCellCentroids''.']);
    end
    dist = sqrt((cy-Gf.cells.centroids(search,1)).^2+(cx-Gf.cells.centroids(search,2)).^2);
    [~,ind] = min(dist);
    centers(i) = search(ind);
end
CG.cells.centers = centers;
end

function pt = geometricMedian(pts)
    pt = mean(pts);
    len = @(v) sqrt(sum(v.^2, 2));
    for i = 1:10
        d = len(bsxfun(@minus, pts, pt));
        pt = sum(bsxfun(@rdivide, pts, d), 1)./sum(1./d, 1);
    end
end
