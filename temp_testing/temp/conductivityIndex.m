function CI = conductivityIndex(G,fracplanes,fraCells,varargin)

opt = struct('tolerance', 0, 'gridDens', 10, 'uniqueTolerance', 1e-2);
opt = merge_options(opt, varargin{:});
gdens = opt.gridDens;
if numel(gdens)==1, gdens = repmat(gdens,3,1); end
pdens = 1/10./gdens;
utol = opt.uniqueTolerance;
Gf = struct;
CI = cell(numel(fracplanes),1);

for i = 1:numel(fracplanes)
    nc = G.nodes.coords(unique(gridCellNodes(G,fraCells{i,1})),:);
    cc = G.cells.centroids(fraCells{i,1},:);

    gridpoints = generatePointsOnPlane(fracplanes(i).points, 'normal', ...
        fracplanes(i).normal, 'ratio',1/gdens);
    tripoints = generatePointsOnPlane(fracplanes(i).points, 'normal', ...
        fracplanes(i).normal, 'ratio',pdens);
    tripoints = uniquetol([gridpoints;tripoints],utol,'ByRows',true);
    
    tri = delaunayTriangulation([cc;nc;tripoints]);
    
    fieldname = ['Plane',num2str(i)];
    Gf(i).(
    
    
end