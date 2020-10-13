function intFun = makeDivQadratureRule(G, cells, degree)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));

nodes = reshape(nodes, 2, [])';

x0 = G.nodes.coords(nodes(:,1),:);
x1 = G.nodes.coords(nodes(:,2),:);

[xr, w] = getQuadratureRule(degree+1, 1);
xr = reshape(xr, 1, 1, []);
w = reshape(w, 1, []);

x = ((x1 - x0).*xr +  (x1 + x0))/2;

ncf = diff(G.cells.facePos);
sign = 1 - 2*(G.faces.neighbors(faces,1) ~= rldecode(cells, ncf(cells), 1));
nx = G.faces.normals(faces,1).*sign/2;

intFun = @(p) integrate(p, cells, faces, w.*nx, x, ncf, G);

end

function val = integrate(p, cells, faces, w, x, ncf, G)

    pz = Polynomial([0,0], 0);
    px = Polynomial([1,0], 1);
    
    q = px*p./(p.k(:,1)+1);
    pp = p./(p.k(:,1)+1);
    
    c = rldecode(cells, ncf(cells), 1);
    
    val = zeros(numel(faces), size(w,2));
    for wNo = 1:size(w,2)
        xhat = (x(:, :, wNo) - G.cells.centroids(c,:))./(G.cells.diameters(c));
        xhat = (x(:, :, wNo) - G.cells.centroids(c,:));%./(G.cells.diameters(c));
%         val(:,wNo) = px(x(:,:,wNo)).*pp(xhat);
        val(:,wNo) = q(xhat);
    end
    
    val = sum(val.*w, 2);
    val = accumarray(rldecode((1:numel(cells))',ncf(cells),1), val);
    
end
