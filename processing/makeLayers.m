function Gl = makeLayers(G,flayers,nl)
%
nfl = numel(flayers);

% Layered matrix grid
Gl = computeGeometry(makeLayeredGrid(G,nl));
if isfield(G,'cartDims')
    Gl.cartDims = [G.cartDims nl];
end
FG = G.FracGrid;
% Layered fracture grid
cstart = Gl.cells.num+1; nstart = Gl.nodes.num+1; fstart = Gl.faces.num+1;
for i = 1:numel(fieldnames(FG))
    tempfl = makeLayeredGrid(FG.(['Frac',num2str(i)]),nfl);
    tempfl.nodes.coords(:,3) = tempfl.nodes.coords(:,3) + min(flayers) - 1;
%     tempfl.cells.centroids(:,3) = tempfl.cells.centroids(:,3) + min(flayers) - 1;
%     tempfl.faces.centroids(:,3) = tempfl.faces.centroids(:,3) + min(flayers) - 1;
    tempfl = computeGeometry(tempfl);
    tempfl.cells.start = cstart;
    tempfl.faces.start = fstart;
    tempfl.nodes.start = nstart;
    Gl.FracGrid.(['Frac',num2str(i)]) = tempfl;
    cstart = tempfl.cells.start + tempfl.cells.num;
    fstart = tempfl.faces.start + tempfl.faces.num;
    nstart = tempfl.nodes.start + tempfl.nodes.num;
end

return