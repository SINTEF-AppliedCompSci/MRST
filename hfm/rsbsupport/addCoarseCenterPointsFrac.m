function CGf = addCoarseCenterPointsFrac(CGf,varargin)
% addCoarseCenterPointsFrac adds coarse nodes to the fracture coarse
% grid CGf. This function can be modified to improve coarse node selection
% inside fractures.
%
% SYNOPSIS:
%   CGf = addCoarseCenterPointsFrac(CGf)
%   CGf = addCoarseCenterPointsFrac(CGf, 'pn1', 'pv1', ...)
%
% REQUIRED PARAMETERS:
%
%   CGf  - Fracture coarse grid (supplied by 'generateCoarseGrid') with
%         geometry information (computed through 'coarsenGeometry').
%
% OPTIONAL PARAMETERS:
%
%   option - option specifies the type of coarse grid points to use for
%            computing the coarse node location. Possible self-explanatory
%            values are: 
%
%            (a) 'useCoarseFaceCentroids'
%            (b) 'useCoarseCellCentroids'
%            (c) 'useFineCellCentroids'
%            (d) 'useCoarseCellEndPoints' - 2D grids only
%               
%            passed as character arrays.
%
%   meantype - type of mean, of the points specified by option, to use for
%              computing the coarse node location. Valid input values
%              include 'geometric' and 'arithmetic' passed as character
%              arrays. This option is not useful if
%              'useCoarseCellEndPoints' is passed as the option
%
% RETURNS:
%   CGf - Fracture coarse grid with CGf.cells.centers added as a new field.
%
% SEE ALSO:
%   storeFractureInteractionRegion

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


Gf = CGf.parent;
opt = struct('option', 'useCoarseFaceCentroids','meantype','geometric'); 
opt = merge_options(opt, varargin{:});

if strcmp(opt.option,'useCoarseCellEndPoints') && CGf.griddim>2
    opt.option = 'useCoarseFaceCentroids';
end

centers = zeros(CGf.cells.num,1);
for i = 1:CGf.cells.num
    search = find(CGf.partition==i);
    if strcmp(opt.option, 'useCoarseFaceCentroids')
        faces = gridCellFaces(CGf,i);
        pts = CGf.faces.centroids(faces,:);
        if strcmp(opt.meantype,'geometric')
            gm = geometricMedian(pts);
        else
            gm = mean(pts);
        end
        C = gm;
    elseif strcmp(opt.option, 'useCoarseCellCentroids')
        C = CGf.cells.centroids(i,:);
    elseif strcmp(opt.option, 'useFineCellCentroids')
        pts = Gf.cells.centroids(search,:);
        if strcmp(opt.meantype,'geometric')
            gm = geometricMedian(pts);
        else
            gm = mean(pts);
        end
        C = gm;
    elseif strcmp(opt.option, 'useCoarseCellEndPoints')
        nbrnum = zeros(numel(search),1);
        for j = 1:numel(search)
            n = getCellNeighbors(Gf, search(j));
            nbrnum(j) = numel(n);
        end
        [~,loc] = min(nbrnum);
        centers(i) = search(loc);
        continue;
    else
        error(['Wrong value for key ''option'' in argument list. Can be either one of ',...
            '''useCoarseFaceCentroids'', ''useCoarseCellCentroids'' or ''useFineCellCentroids''.']);
    end
    diff = repmat(C,numel(search),1)-Gf.cells.centroids(search,:);
    dist = sqrt(sum(diff.^2,2));
    [~,ind] = min(dist);
    centers(i) = search(ind);
end
CGf.cells.centers = centers;
return

function pt = geometricMedian(pts) % Olav
    pt = mean(pts);
    len = @(v) sqrt(sum(v.^2, 2));
    for i = 1:10
        d = len(bsxfun(@minus, pts, pt));
        pt = sum(bsxfun(@rdivide, pts, d), 1)./sum(1./d, 1);
    end
return