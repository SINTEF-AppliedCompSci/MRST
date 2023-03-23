function tests = mpsaTest
% Unit tests for MPSA
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
Nx = [2 2];  Nd = numel(Nx);
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

function translationTest(testCase)
mu = 100;
lambda = 100;
phi = 0;
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
    
    constit = shear_normal_stress(Nc, Nd, mu*ones(Nc,1), lambda*ones(Nc,1), phi*ones(Nc,1));
    md = mpsa(G,constit,[],'bc',bc);
    
    % Translation
    r = rand(1,Nd);
    db = ones(Nf,1) * r;
    xan = ones(Nc,1) * r;
    boundVal = -md.div * md.boundStress * reshape(db',[],1);
    
    x = reshape((md.A \ boundVal)',Nd,[])';
    assert(max(max(abs(x - xan)))<sqrt(eps),'MPSA failed on translation');
    
    stress = md.stress * reshape(x',[],1) + md.boundStress * reshape(db',[],1);
    
    assert(norm(stress) < sqrt(eps),'MPSA gave stress for translation')
end
end


function rotationTest(testCase)
mu = 100;
lambda = 100;
phi = 0;
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
    xn = G.nodes.coords;
    
    if Nd == 2
        ang = 180*rand(1);
        rot = @(x) [x(:,1) * cosd(ang) - x(:,2) * sind(ang), x(:,1) * sind(ang) + x(:,2) * cosd(ang)];
    else
        ang = 180 * rand(3,1);
        rotx = @(x) [x(:,1) , x(:,2) * cosd(ang(1)) - x(:,3) * sind(ang(1)), x(:,2) * sind(ang(1)) + x(:,3) * cosd(ang(1))];
        roty = @(x) [x(:,1) * cosd(ang(2)) - x(:,3) * sind(ang(2)), x(:,2) ,  x(:,1) * sind(ang(2)) + x(:,3) * cosd(ang(2))];
        rotz = @(x) [x(:,1) * cosd(ang(3)) - x(:,2) * sind(ang(3)), x(:,1) * sind(ang(3)) + x(:,3) * cosd(ang(3)), x(:,3)];
        rot = @(x) rotz(roty(rotx(x)));
    end
    
    db  = rot(xf) - xf;
    xan = rot(xc) - xc;
   
    constit = shear_normal_stress(Nc, Nd, mu*ones(Nc,1), lambda*ones(Nc,1), phi*ones(Nc,1));
    md = mpsa(G,constit,[],'bc',bc,'eta',0);
    x = reshape(-(md.A) \ md.div * md.boundStress * reshape(db',[],1),Nd,[])';
    assert(max(max(abs(x - xan)))<sqrt(eps),'MPSA failed on rotation');
end
end


function uniformStrainTest(testCase)
mu = 100;
lambda = 100;
phi = 0;
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
    xn = G.nodes.coords;

    % Uniform stress
    r = rand(Nd,1);
    db = bsxfun(@times,xf,r');
    xan = bsxfun(@times,xc,r');
   
  
    constit = shear_normal_stress(Nc, Nd, mu*ones(Nc,1), lambda*ones(Nc,1), phi*ones(Nc,1));
    md = mpsa(G,constit,[], 'bc',bc);
     
        x = reshape(-full(md.A) \ md.div * md.boundStress * reshape(db',[],1),Nd,[])';

    assert(max(max(abs(x - xan)))<sqrt(eps),'MPSA failed on uniform stress');
    stress = md.stress * reshape(x',[],1) + md.boundStress * reshape(db',[],1);
    
    san = 2 * mu * r' + lambda * sum(r);
    san = (ones(Nf,1) * san).*G.faces.normals;
    s = reshape(stress',Nd,[])';
    
    assert(max(max(abs(s - san))) < sqrt(eps),'MPSA failed on stress for uniform stretching')
    
end
end