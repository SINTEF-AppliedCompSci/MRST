function vc = cellAverage(G, vf)
% compute cell-average vc of field vf as average of adjecent faces
if numel(vf) ~= G.faces.num % assume interior
    ix   = ~any(G.faces.neighbors==0, 2); %index to interor faces
    assert(numel(vf)==nnz(ix))
    v    = zeros(G.faces.num,1); 
    v(ix) = vf;
else
    ix = true(G.faces.num,1);
    v  = vf;
end
cellno = gridCellNo(G); 
numCF   = accumarray(cellno, ix(G.cells.faces(:,1))); % number of faces per cell
sumVals = accumarray(cellno, v(G.cells.faces(:,1))); % sum of trans per cell
vc = sumVals./numCF;
end