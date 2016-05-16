function [fraCells, remove] = markcells(G, fracplanes, varargin)
% markcells returns the indices of matrix cells that are connected with
% each fracture.

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
opt = struct('inPolygonTolerance',     0, ...
             'Sign'              ,     [], ...
             'GlobTri'           ,     struct('Tri', [], 'map', []));
         
opt = merge_options(opt, varargin{:});
Gtri = opt.GlobTri.Tri;
map = opt.GlobTri.map;
maxd = 0.5*max(G.cells.volumes).^(1/3);
if isempty(opt.Sign)
    Sign = establishSign(G, fracplanes);
else
    Sign = opt.Sign;
end
remove = [];
fraCells = cell(numel(fracplanes),1);
%
[fnodes,fnmap] = gridFaceNodes(G,1:G.faces.num);
nodesPerFace = diff(fnmap);
%
possible_ratios = 0:0.0001:1;
possible_ratios = possible_ratios(mod(1,possible_ratios)==0);
%
xc = G.cells.centroids(:,1); yc = G.cells.centroids(:,2); zc = G.cells.centroids(:,3);
%
if any(diff(nodesPerFace)) % Hybrid/Composite Grid
    faceNodes = zeros(G.faces.num,max(nodesPerFace));
    for i = 1:G.faces.num
        temp = fnodes(fnmap(i):fnmap(i+1)-1).';
        faceNodes(i,1:numel(temp)) = temp;
    end
else
    faceNodes = reshape(fnodes,nodesPerFace(1),[]).';
end
count = 0; flag = 0;
for i = 1:numel(fracplanes)
    count = count+1;
    %
    ratio = max(G.faces.areas)/(polyArea3D(fracplanes(i).points));
    [~,loc] = min(abs(possible_ratios-ratio));
    if loc>=numel(possible_ratios)-2, loc = loc-4; end
    ratio = possible_ratios(loc+1);
    %
    bsum = Sign(i).NodeSign(:,2:end)>=0; % Points lying on boundary planes
    onplane = Sign(i).NodeSign(:,1) == 0; % Points on plane
    nodesOfInterest  = find(sum([onplane,bsum],2) == size(fracplanes(i).points,1)+1);
    %
    isFace = ismember(faceNodes,nodesOfInterest);
    dfm_cells = unique(G.faces.neighbors(sum(isFace,2)==sum(faceNodes~=0,2),:));
    %
    [p,in,on] = generatePointsOnPlane(fracplanes(i).points, 'ratio', ratio, ...
                'tolerance',opt.inPolygonTolerance);
    %
    p = p(in~=on,:);
    if ~isempty(dfm_cells)
        unitNormal = fracplanes(i).normal/norm(fracplanes(i).normal);
        tdist = fracplanes(i).aperture;
        pminus = p - repmat(unitNormal,size(p,1),1)*tdist;
        pplus = p + repmat(unitNormal,size(p,1),1)*tdist;
        p = [pminus;pplus];
    end
    if ~isempty(Gtri)
        ploc = pointLocation(Gtri,p);
        ploc = ploc(~isnan(ploc));
        allc = unique([map(ploc);dfm_cells]);
        allc = allc(allc~=0);
        xct = xc(allc); yct = yc(allc); zct = zc(allc);
        dd = zeros(size(xct,1),1);
        for j = 1:numel(allc)
            dd(j) = min(sqrt((p(:,3)-repmat(zct(j),size(p,1),1)).^2 + (p(:,2)-repmat(yct(j),size(p,1),1)).^2 + (p(:,1)-repmat(xct(j),size(p,1),1)).^2));
        end
        allc = setdiff(allc,allc(dd>maxd));
        if numel(allc)<=3, flag = 1; end
    end
    if isempty(Gtri) || flag == 1
        flag = 0;
        dd = zeros(G.cells.num,1);
        for j = 1:G.cells.num
            dd(j) = min(sqrt((p(:,3)-repmat(zc(j),size(p,1),1)).^2 + (p(:,2)-repmat(yc(j),size(p,1),1)).^2 + (p(:,1)-repmat(xc(j),size(p,1),1)).^2));
        end
        candidates = find(dd<maxd);
        dd = dd(dd<maxd);
        [~,sortInd] = sort(dd);
        [cells, ~] = getEnclosingCellsByFace(G, p, 'candidates', candidates(sortInd));
        allc = unique([cells;dfm_cells]);
    end
    if numel(allc)<=3
        warning(['Fracture ',num2str(i),' is not a long fracture and ',...
            'will be removed from further calculations.']);
        remove = [remove;i]; %#ok
        count = count-1;
        continue
    end
    fraCells{count,1} = allc;
end

return


    
            