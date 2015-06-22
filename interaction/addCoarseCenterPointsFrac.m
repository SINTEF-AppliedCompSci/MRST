function CGf = addCoarseCenterPointsFrac(CGf,varargin)
% addCoarseCenterPointsFrac(CG) adds coarse nodes to the fracture coarse
% grid CG
%
% OPTIONAL PARAMETERS (supplied in 'key'/'value' pairs ('pn'/'pv' ...)):
%
%   option - option specifies the type of coarse grid points to use for
%            computing the coarse node location. Possible values are:
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

Gf = CGf.parent;
opt = struct('option', 'useCoarseFaceCentroids','meantype','geometric'); 

opt = merge_options(opt, varargin{:});
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