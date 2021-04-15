function I = polyhedronInt(G,cells,f, k)
%Integrates the function f over each cell in cells of grid G, using a
%quadrature rule of precision k.
%
%   SYNOPSIS:
%       I = polygonInt(G, cells, f, k)
%
%   DESCRIPTION:
%       Approximates the integrals
%
%           \int_K f \dx
%
%       over specified cells K of G of using a quadrature rule of
%       precission k. Each cell divided into terahedra, and a map F from
%       reference tetrahedron Tr with vertices V constructed. Using that
%
%           \int_T f \dx = |\det(F)|\int_Tr f(F(y)) \dy,
%
%       the integral can be approximated by the quadrature rule.
%
%   REQUIRED PARAMETERS:
%       G       - MRST grid.
%       cells   - Cells over which to integrate f.
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

[Xq, w, V, vol] = tetrahedronQuadRule(k);
    
Vdiff = V(1:end-1,:) - V(2:end,:);
   
nK = numel(cells);
nq = size(Xq,1);
      
I = zeros(nK,size(f([0,0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);
    
    X = X(reshape(tri',[],1),:);
    ii = mcolon(1:4:4*nTri, 3:4:4*nTri);
    Xdiff = X(ii,:) - X(ii+1,:);

    ii = repmat((1:3)',3*nTri,1);
    add = repmat((0:1:nTri-1)*3,9,1);
    ii = ii + add(:);
    jj = repmat((1:3), 3, 1);
    jj = repmat(jj(:), nTri,1);
    jj = jj(:) + add(:);
    Vd = sparse(ii,jj,repmat(Vdiff(:),nTri,1), 3*nTri, 3*nTri);
    Xdiff = Xdiff';
    Xd = sparse(jj,ii,Xdiff(:), 3*nTri, 3*nTri);
    
    A = Vd\Xd;

    Xb = reshape(X(1:4:end,:)',1,[]);

    b = Xb - repmat(V(1,:),1,nTri)*A;

    [n,~] = size(A);
    
    main = ones(n,1);
    sub1 = main;
    sub1(3:3:end) = 0;
    sub2 = zeros(n,1);
    sub2(1:3:end) = 1;
    Id = spdiags([sub2, sub1, main,sub1(end:-1:1),sub2(end:-1:1)], -2:1:2, n,n);

    Ad = full(A(Id == 1));
    Ad = reshape(Ad,3,3,[]);
    D = squeeze(Ad(1,1,:).*Ad(2,2,:).*Ad(3,3,:) + ...
                Ad(2,1,:).*Ad(3,2,:).*Ad(1,3,:) + ...
                Ad(3,1,:).*Ad(1,2,:).*Ad(2,3,:) - ...
                Ad(3,1,:).*Ad(2,2,:).*Ad(1,3,:) - ...
                Ad(2,1,:).*Ad(1,2,:).*Ad(3,3,:) - ...
                Ad(1,1,:).*Ad(3,2,:).*Ad(2,3,:));
    D = abs(D);
    
    XhatTmp = bsxfun(@plus, repmat(Xq, 1, nTri)*A, b);
    Xhat = zeros(nq*nTri,3);
    Xhat(:,1) = reshape(XhatTmp(:,1:3:3*nTri),[],1);
    Xhat(:,2) = reshape(XhatTmp(:,2:3:3*nTri),[],1);
    Xhat(:,3) = reshape(XhatTmp(:,3:3:3*nTri),[],1);

    I(i,:) = vol*repmat(w,1,nTri).*rldecode(D,nq*ones(nTri,1),1)'*f(Xhat);

end

end