function G = fracMatrixConnections(G,Gfrac,CItot,possible_cells,area,varargin)
% fracMatrixConnections assigns a "non-neighboring connection (NNC)" status
% to each fracture-matrix connection and also assigns a transmissibility to
% each NNC given the total conductivity index of the fracture. See Lee et
% al, Water Resources Research, 2001 or SPE-65095-PA, Lee et al, 2000.

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

opt = struct('GlobTri',struct('Tri',[],'map',[]));
         
opt = merge_options(opt, varargin{:});
Gtri = opt.GlobTri.Tri;
map = opt.GlobTri.map;

vr = Gfrac.cells.volumes./sum(Gfrac.cells.volumes);
[cn,cpos] = gridCellNodes(Gfrac,1:Gfrac.cells.num);
if isempty(Gtri)
    cells = getEnclosingCellsByFace(G,Gfrac.nodes.coords);
else
    try
        cells = map(pointLocation(Gtri,Gfrac.nodes.coords));
    catch
        cells = getEnclosingCellsByFace(G,Gfrac.nodes.coords);
    end
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
    