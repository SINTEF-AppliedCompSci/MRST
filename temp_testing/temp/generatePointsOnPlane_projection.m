function copl = generatePointsOnPlane_projection(G,fracplanes,fraCells,varargin)

opt = struct('tolerance', 0, 'linspaceParts', 100, 'ratio', 1e-2);
opt = merge_options(opt, varargin{:});
pdens = opt.linspaceParts;
if numel(pdens)==1, pdens = repmat(pdens,3,1); end


for i = 1:numel(fracplanes)
    nc = G.nodes.coords(unique(gridCellNodes(G,fraCells{i,1})),:);
    ma = max(nc);
    mi = min(nc);
    xx = linspace(mi(1),ma(1),pdens(1));
    yy = linspace(mi(2),ma(2),pdens(2));
    zz = linspace(mi(3),ma(3),pdens(3));
    [X,Y,Z] = meshgrid(xx,yy,zz);
    pbox = [X(:),Y(:),Z(:)];
    pbox = [pbox;pbox(1,:)];
    
    points = fracplanes(i).points;
    normal = fracplanes(i).normal;
    
    A = normal(1); B = normal(2); C = normal(3);
    D = -dot(normal,points(1,:));

    
    basisCoefficients= bsxfun(@minus,pbox,points(1,:))*null(normal);
    copl = bsxfun(@plus,basisCoefficients*null(normal).', points(1,:));
    copl = uniquetol(copl,1e-2,'ByRows',true);
    
end