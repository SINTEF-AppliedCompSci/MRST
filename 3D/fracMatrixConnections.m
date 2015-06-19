function G = fracMatrixConnections(G,Gfrac,CItot,possible_cells,area,varargin)

opt = struct('GlobTri',struct('Tri',[],'map',[]));
         
opt = merge_options(opt, varargin{:});
Gtri = opt.GlobTri.Tri;
map = opt.GlobTri.map;

vr = Gfrac.cells.volumes./sum(Gfrac.cells.volumes);
[cn,cpos] = gridCellNodes(Gfrac,1:Gfrac.cells.num);
if isempty(Gtri)
    cells = getEnclosingCellsByFace(G,Gfrac.nodes.coords);
else
    cells = map(pointLocation(Gtri,Gfrac.nodes.coords));
end
for i = 1:Gfrac.cells.num
    nodes = cn(cpos(i):cpos(i+1)-1);
    connCells = unique(cells(nodes));
    connCells = connCells(ismember(connCells,possible_cells));
    G.nnc.cells = [G.nnc.cells; connCells, ...
                   repmat(i+Gfrac.cells.start-1,numel(connCells),1)];
    G.nnc.CI = [G.nnc.CI; repmat(CItot*vr(i)/numel(connCells),numel(connCells),1)];
    G.nnc.area = [G.nnc.area; repmat(area*vr(i)/numel(connCells),numel(connCells),1)];
end
return
    