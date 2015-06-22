function [fraCells, remove] = markcells_n(G, fracplanes, varargin)
opt = struct('pointDens',100);
         
opt = merge_options(opt, varargin{:});
pdens = 1/opt.pointDens;

remove = [];
fraCells = cell(numel(fracplanes),1);
nc = G.nodes.coords;
cc = G.cells.centroids;
numc = G.cells.num;
count = 0;
for i = 1:numel(fracplanes)
    count = count+1;
    tripoints = generatePointsOnPlane(fracplanes(i).points, 'normal', ...
        fracplanes(i).normal, 'ratio', pdens);
    tri = delaunayTriangulation([cc;tripoints;nc]); % 1:G.cells.num entries will be centroids
    poi = [1:numc,numc+1:numc+size(tripoints,1)]';
    T = tri.ConnectivityList;
    check = ismember(T,poi);
    T = T(sum(check,2)==4,:);
    T = unique(T(T<=numc));
    if numel(T)<=3
        warning(['Fracture ',num2str(i),' is not a long fracture and ',...
            'will be removed from further calculations.']);
        remove = [remove;i]; %#ok
        count = count-1;
        continue
    end
    fraCells{count,1} = T;
end