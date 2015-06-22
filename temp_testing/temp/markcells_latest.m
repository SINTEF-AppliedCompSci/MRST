function [fraCells, remove] = markcells(G, GlobTri, fracplanes, varargin)

opt = struct('inPolygonTolerance',       0, ...
             'Sign'              ,      [], ...
             'distTolerance'     ,    eps*5 );
         
opt = merge_options(opt, varargin{:});

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
Gtri = GlobTri.Tri;
map = GlobTri.map;
%
possible_ratios = 0:0.0001:1;
possible_ratios = possible_ratios(mod(1,possible_ratios)==0);

if any(diff(nodesPerFace)) % Hybrid/Composite Grid
    faceNodes = zeros(G.faces.num,max(nodesPerFace));
    for i = 1:G.faces.num
        temp = fnodes(fnmap(i):fnmap(i+1)-1).';
        faceNodes(i,1:numel(temp)) = temp;
    end
else
    faceNodes = reshape(fnodes,nodesPerFace(1),[]).';
end
count = 0;
for i = 1:numel(fracplanes)
    count = count+1;
    %
    ratio = max(G.faces.areas)/(10*polyArea3D(fracplanes(i).points));
    [~,loc] = min(abs(possible_ratios-ratio));
    ratio = possible_ratios(loc);
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
%     [~,disit] = iscoplanar(p,'normal',fracplanes(i).normal,...
%         'distTolerance', opt.distTolerance, 'point', fracplanes(i).points(1,:));
    %
    p = p(in~=on,:); % & disit,:);
    if ~isempty(dfm_cells)
        unitNormal = fracplanes(i).normal/norm(fracplanes(i).normal);
        tdist = fracplanes(i).aperture/2;
        pminus = p - repmat(unitNormal,size(p,1),1)*tdist;
        pplus = p + repmat(unitNormal,size(p,1),1)*tdist;
        p = [pminus;pplus];
    end
%     ploc = pointLocation(Gtri,p);
%     ploc = ploc(~isnan(ploc));
%     cells = unique(map(ploc));
    [cells, ~] = getEnclosingCellsByFace(G, p);
    allc = unique([cells;dfm_cells]);
    allc = allc(allc~=0);
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


    
            