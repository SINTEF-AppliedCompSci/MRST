function I = polygonInt(G, cells, f, k)
%Integrates the function f over each cell in cells of grid G, using a
%quadrature rule of precision k.
%
%   SYNOPSIS:
%       I = polygonInt(G, cells, f, k)
%
%   DESCRIPTION:
%       Approximates the integrals
%           \int_K f \dx
%       over specified cells K of G of using a quadrature rule of
%       precission k. Each cell is trangulated, and a map F from reference
%       triangle with vertices (0,0), (1,0) and (0,1) is constructed. Using
%       that
%
%           \int_K f \dx = |\det(F)|\int_T f(F(y)) \dy,
%
%       the integral can be approximated by the quadrature rule.
%
%   REQUIRED PARAMETERS:
%       G       - MRST grid.
%       cells   - Cells over which to integrate f.
%       f       - Integrand.
%       k       - Precission of quadrature rule.
%
%   RETURNS:
%       I       - Approximated solution to the integral.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.


[Xq, w, ~, vol] = triangleQuadRule(k);

nq = size(Xq,1);

nK = numel(cells);

I = zeros(nK,size(f([0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    
    X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);

    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    Ad = reshape(A',2,2,[]);
    D = abs(squeeze(Ad(1,1,:).*Ad(2,2,:) - Ad(1,2,:).*Ad(2,1,:)));
    
    ii = repmat((1:2*nTri)',2,1);
    jj = 1:2*nTri;
    jj = reshape(jj,2,[])';
    jj = repmat(jj(:)',2,1);
    jj = jj(:);
    
    A = sparse(ii, jj, A(:), 2*nTri, 2*nTri);
    
    XhatTmp = repmat(Xq,1,nTri)*A + repmat(reshape(bA',1,[]),nq,1);
    Xhat = zeros(nq*nTri,2);
    Xhat(:,1) = reshape(XhatTmp(:,1:2:end),[],1);
    Xhat(:,2) = reshape(XhatTmp(:,2:2:end),[],1);
    
    I(i,:) = vol*repmat(w,1,nTri).*(rldecode(D,nq*ones(nTri,1),1))'*f(Xhat);
     
end

end