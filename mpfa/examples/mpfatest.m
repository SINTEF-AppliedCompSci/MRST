%% Simple test where we compare the different implementation: legacy, tensor assembly, with and without iteration over blocks (necessary for large systems)
%
% The test is aimed at Neumann boundary conditions (no flow) and well injection.
%

clear all

mrstModule add mimetic mpfa incomp

isverbose = true;
eta       = 1/3;
blocksize = 10;

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 10; ny = 10;
G = cartGrid([nx, ny]);
G = twister(G, 0.1);
G = computeGeometry(G);
nc = G.cells.num;

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
% We make a non diagonal rock tensor
rock = makeRock(G, 1e-3*darcy, 1);

% injcells  =  1 : 3;
% prodcells =  nc : -1 : (nc - 3);
injcells  =  1 : 3;
prodcells =  nc : -1 : (nc - 2);
radius    = 0.1;

rate = 1e-4*meter^3/day;
bhp  = 1*barsa;

W = addWell([], G, rock, injcells, 'Comp_i', 1, 'Type', 'rate', 'Val', rate, ...
            'sign', 1, 'Radius', radius);
W = addWell(W, G, rock, prodcells, 'Comp_i', 1, 'Type', 'bhp', 'Val', bhp, ...
            'sign', -1, 'Radius', radius);

fluid = initSingleFluid('mu' , 1 , ...
                        'rho', 1*kilogram/meter^3);

gravity off

clear mpfastructs states pressures fluxes titles
caseno = 1;

state0 = initResSol(G, 0, 1);

%% mpfa - legacy
tic
T_mpfa = computeMultiPointTrans(G, rock, 'eta', eta);
texec = toc;
states{caseno}    = incompMPFA(state0, G, T_mpfa, fluid, 'wells', W);
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno}    = states{caseno}.flux;
titles{caseno}    = 'mpfa - legacy';
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
caseno = caseno + 1;

%% mpfa - tensor assembly
tic
opts = {'eta'              , eta      , ...
        'verbose'          , isverbose, ...
        'useTensorAssembly', true     , ...
        'neumann'          , true};
mpfastructs{caseno} = computeMultiPointTrans(G, rock, opts{:});
texec = toc;
% Options for incompMPFA for tensor assembly call
incompopts = {'useTensorAssembly', true, ...
              'wells'            , W   , ...
              'outputFlux'       , true};
states{caseno}    = incompMPFA(state0, G, mpfastructs{caseno}, [], incompopts{:});
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno}    = states{caseno}.flux;
titles{caseno}    = 'mpfa - Neumann';
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
caseno = caseno + 1;

%% mpfa - tensor assembly with iterations on blocks (necessary for large system)
tic
opts = {'eta'              , eta      , ...
        'verbose'          , isverbose, ...
        'useTensorAssembly', true     , ...
        'blocksize'        , blocksize, ...
        'neumann'          , true};
mpfastructs{caseno} = computeMultiPointTrans(G, rock, opts{:});
texec = toc;
states{caseno} = incompMPFA(state0, G, mpfastructs{caseno}, [], incompopts{:});
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno} = states{caseno}.flux;
titles{caseno} = 'mpfa - Neumann - block';
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
caseno = caseno + 1;


%% Plotting
for i = 1 : numel(pressures)
    figure(i)
    clf
    plotCellData(G, pressures{i} ./ barsa());
    title(titles{i});
    colorbar
end


%% check flux computations by computing mass directly mass conservation in
%% each cell

nc = G.cells.num;
nf = G.faces.num;

N = G.faces.neighbors;

cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
cellfacetbl.faces = G.cells.faces(:, 1);
cellfacetbl = IndexArray(cellfacetbl);

celltbl.cells = (1 : nc)';
celltbl = IndexArray(celltbl);

facetbl.faces = (1 : nf)';
facetbl = IndexArray(facetbl);

fno = cellfacetbl.get('faces');
cno = cellfacetbl.get('cells');
sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;

prod = TensorProd();
prod.tbl1 = cellfacetbl;
prod.tbl2 = facetbl;
prod.tbl3 = celltbl;
prod.reducefds = {'faces'};
prod = prod.setup();

div_T = SparseTensor();
div_T = div_T.setFromTensorProd(sgn, prod);
div = div_T.getMatrix();

clear qs;
for i = 1 : numel(fluxes)
    qs{i} = div*(fluxes{i});
    figure
    plotCellData(G, qs{i});
    title(['mass conservation - ', titles{i}]);
    colorbar
    figure
    plot(qs{i}, '*');
    title(['mass conservation - ', titles{i}]);
    xlabel('cell number');
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
