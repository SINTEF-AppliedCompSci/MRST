mrstModule add dg upr

%%

usePebi = true;
n       = 20;

if usePebi
    % Create 2D PEBI grid
    bdr = [-1, -1;
            1, -1;
            1,  1;
            -1, 1];
    G2D = pebiGrid(1/sqrt(n), [2,2], 'polybdr', bdr);
    G2D = computeVEMGeometry(G2D);
    G2D = computeCellDimensions(G2D);
    G2D.name = '2D PEBI';
    % Create 3D PEBI grid
    bdr = [-1, -1, -1;
            1, -1, -1;
            1,  1, -1;
           -1,  1, -1;
           -1, -1,  1;
            1, -1,  1;
            1,  1,  1;
           -1,  1,  1];
    bdr = delaunayTriangulation(bdr);
    pts = rand(100, 3)*2 - 1;
    G3D = clippedPebi3D(pts, bdr);
    G3D = computeVEMGeometry(G3D);
    G3D = computeCellDimensions(G3D);
    G3D.name = '3D PEBI';
else
    % Create 2D Cartesian grid
    G2D = cartGrid([n,n], [2,2]);
    G2D.nodes.coords = G2D.nodes.coords - 1;
    G2D = computeVEMGeometry(G2D);
    G2D = computeCellDimensions(G2D);
    G2D.name = '2D Cart';
    % Create 3D Cartesian grid
    G3D = cartGrid([n,n,n], [2,2,2]);
    G3D.nodes.coords = G3D.nodes.coords - 1;
    G3D = computeVEMGeometry(G3D);
    G3D = computeCellDimensions(G3D);
    G3D.name = '3D Cart';
end
grids = {G2D, G3D};

%%

tol  = 1e-12;
kMax = 6;
px = {Polynomial([1,0], 1), Polynomial([1,0,0], 1)};
py = {Polynomial([0,1], 1), Polynomial([0,1,0], 1)};
pz = {Polynomial([0,0], 0), Polynomial([0,0,1], 1)};

%%


for gNo = 1:numel(grids)

    G = grids{gNo};
    fprintf(['\nTesting cubatures on ', G.name, '...\n']);
    
    for k = 1:kMax
        
        if k == 1
            stndrdth = 'st';
        elseif k == 2
            stndrdth = 'nd';
        elseif k == 3
            stndrdth = 'rd';
        else
            stndrdth = 'th';
        end
    
        if G.griddim == 2
            cubVol  = TriangleCubature(G, k, []);
            cubVol = MomentFitting2DCubature(G, k, []);
            cubSurf = LineCubature(G, k, []);
        else
            cubVol  = TetrahedronCubature(G, k, []);
            cubSurf = TriangleCubature(G, k, []);
        end

        basis  = dgBasis(G.griddim, k, 'legendre');
        nDof   = basis.nDof;
        sol    = zeros(nDof,1);
        sol(1) = sum(G.cells.volumes);

        % Test volume cubatures
        [W, x] = cubVol.getCubature((1:G.cells.num)', 'volume');    
        I = zeros(nDof,1);
        for dofNo = 1:nDof
            I(dofNo) = sum(W*basis.psi{dofNo}(x));
        end
        assert(all(abs(I - sol) < tol), ...
            [num2str(k), stndrdth, ' order ', class(cubVol), ...
                                                   ' failed on ', G.name]);
        fprintf([num2str(k), stndrdth, ' order ', class(cubVol), ...
                                        ' successfull on ', G.name, '\n']);
       
        % Test surface cubatures
        [W, x, w, ~, fNo] = cubSurf.getCubature((1:G.cells.num)', ...
                                         'surface', 'outwardNormal', true);
        I = zeros(nDof,1);
        n = G.faces.normals./G.faces.areas;
        n = n(fNo,:);
        for dofNo = 1:nDof
            psi = basis.psi{dofNo};
            % Use divergence theorem to convert to surface integral
            psix = cell(G.griddim,1);
            for dNo = 1:G.griddim
                p = psi;
                p.k(:,dNo) = p.k(:,dNo) + 1;
                p.w = p.w./p.k(:,dNo);
                psix{dNo} = p;
            end
            psix = [psix{:}];
            
            f        = sum(psix(x).*n,2)/G.griddim;
            I(dofNo) = sum(W*f);
        end
        stth = 'th';
        if k == 1; stth = 'st'; end
        assert(all(abs(I - sol) < tol), ...
            [num2str(k), stndrdth, ' order ', class(cubSurf), ...
                                                   ' failed on ', G.name]);
        fprintf([num2str(k), stndrdth, ' order ', class(cubSurf), ...
                                        ' successfull on ', G.name, '\n']);
        
    end
    
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
