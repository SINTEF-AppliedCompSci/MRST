function tests = mpfaTest
% Unit tests for MPFA.
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
grids = {};

% 2D Cartesian, no perturbations
Nx = [4 3];  Nd = numel(Nx);
g = cartGrid(Nx); Nn = g.nodes.num;
grids{1} = computeGeometry(g);

% Perturbations
g.nodes.coords = g.nodes.coords + 0.5 * rand(Nn,Nd);
grids{2} = computeGeometry(g);

% 3D Cartesian, no perturbations
Nx = [3 3 3];  Nd = numel(Nx);
g = cartGrid(Nx); Nn = g.nodes.num;
grids{3} = computeGeometry(g);

% Perturbations
g.nodes.coords = g.nodes.coords + 0.5 * rand(Nn,Nd);
grids{4} = computeGeometry(g);

% Triangular grid
Nx = [5 5]; Nd = numel(Nx);
[X,Y] = meshgrid(linspace(0,Nx(1),Nx(1)+1), linspace(0,Nx(2),Nx(2)+1));
p     = [X(:), Y(:)];
t     = delaunayn(p);
g     = triangleGrid(p, t); Nn = g.nodes.num;
grids{5} = computeGeometry(g);

% Perturbed triangles
g.nodes.coords = g.nodes.coords + 0.5 * rand(Nn,Nd);
grids{6} = computeGeometry(g);

% Tetrahedral grid
Nx = [2 2 2]; Nd = numel(Nx);
[X,Y,Z] = meshgrid(linspace(0,Nx(1),Nx(1)+1), linspace(0,Nx(2),Nx(2)+1),linspace(0,Nx(3),Nx(3)+1));
p     = [X(:), Y(:), Z(:)];
t     = delaunayn(p);
g     = tetrahedralGrid(p, t); Nn = g.nodes.num;
grids{7} = computeGeometry(g);

% Perturbed triangles
g.nodes.coords = g.nodes.coords + 0.5 * rand(Nn,Nd);
grids{8} = computeGeometry(g);

testCase.TestData.grids = grids;

end


function dirichletBoundaryTest(testCase)

grids = testCase.TestData.grids;

for iter1 = 1 : numel(grids);
    % Not the most beautiful of setups, we could probably have used setup somehow
    G = grids{iter1};
    Nc = G.cells.num;
    Nf = G.faces.num;
    Nd = G.griddim;
    
    xf = G.faces.centroids;
    xc = G.cells.centroids;
    
    % Pressure boundary conditions;
    bc = addBC([],find(any(G.faces.neighbors == 0,2)),'pressure',0);
    
    fd = mpfa(G,struct('perm',ones(Nc,1)),[],'bc',bc,'invertBlocks','matlab');
    
    % First test no flow
    ub = ones(Nf,1);
    xan = ones(Nc,1);
    x = -fd.A \ (fd.div * fd.boundFlux * ub);
    assert(norm(x - xan) < sqrt(eps),'MPFA on no flow failed')
    
    % Then test uniform flow in x-direction
    ub = xf(:,1);
    xan = xc(:,1);
    x = -fd.A\fd.div * fd.boundFlux * ub;
    
    assert(norm(x - xan) < sqrt(eps),'MPFA on no flow in x-direction failed')
    
    % Flow in x+y direction
    r = rand(Nd,1);
    ub = xf * r;
    xan = xc * r;
    x = -fd.A\fd.div * fd.boundFlux * ub;
    assert(norm(x - xan) < sqrt(eps),'MPFA on no flow in x + y-direction failed')
    
    flux =  fd.F * x + fd.boundFlux * ub;
    fan = sum(bsxfun(@times,G.faces.normals,-r'),2);
    
    assert(norm(flux - fan) < sqrt(eps),'MPFA flux error on random field')
    
    
end
end

function mixedBoundaryTest(testCase)
grids = testCase.TestData.grids;

for iter1 = 7 : numel(grids);
    % Not the most beautiful of setups, we could probably have used setup somehow
    G = grids{iter1};
    Nc = G.cells.num;
    Nf = G.faces.num;
    
    xf = G.faces.centroids;
    xc = G.cells.centroids;
    
     if   any(strcmpi(G.type,'tensorGrid'))
        
        bc = fluxside([],G,'xmin',1);
        bc = pside(bc,G,'xmax',1);
        bc = pside(bc,G,'ymin',1);
        bc = pside(bc,G,'ymax',1);
        
        fd = mpfa(G,struct('perm',ones(Nc,1)),[],'bc',bc);
        
        fb = ones(Nf,1);
        
        x = -fd.A \ (fd.div * fd.boundFlux * fb);
        a=[];
     end
end

end
