%% Linear pressure test
%
%
% The MPFA method is exact for linear pressure field. We impose boundary
% condition such that the exact solution is linear. The grid is a twisted grid
% made from a Cartesian grid.
%
% We check the three implementations (Legacy, standard, block assembly)


mrstModule add ad-core ad-props incomp mrst-gui mpfa

clear all

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
dimcase = 2;
switch dimcase
  case 2
    nx = 10; ny = 10;
    G = cartGrid([nx, ny]);
  case 3
    nx = 4; ny = 3; nz = 5;
    G = cartGrid([nx, ny, nz]);
end
G = twister(G, 0.1); 
G = computeGeometry(G);
nc = G.cells.num;

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
% We make a non diagonal rock tensor
% rock = makeRock(G, 1e-3*darcy, 1);
rock = makeRock(G, 1, 1);

fluid = initSingleFluid('mu' , 1, ...
                        'rho', 1);

gravity off

%% setup bc for pressure
dim = G.griddim;
zface = G.faces.centroids(:, dim);
maxz = max(zface);
minz = min(zface);
bottomfaces = find(G.faces.centroids(:, dim) > (maxz - 1e-4));
topfaces = find(G.faces.centroids(:, dim) < (minz + 1e-4));

nbf = numel(bottomfaces);
ntf = numel(topfaces);
bc = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
bc.face = [bottomfaces; topfaces];
bc.type = repmat({'pressure'}, [1, numel(bc.face)]);
value = zeros(nbf + ntf, 1);
value(1 : nbf) = 1*barsa;
bc.value = value;


%% Pressure run
titles = {};
z = G.cells.centroids(:, dim);
eta = 1/3;
blocksize = 10;

state0 = initResSol(G, 0, 1);

clear vecs fluxes
caseno = 1;

%% mpfa - Legacy implementation

T_mpfa = computeMultiPointTrans(G, rock, 'eta', eta);
state = incompMPFA(state0, G, T_mpfa, fluid, 'bc', bc);
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - legacy';
caseno         = caseno + 1;

%% mpfa - Tensor Assembly implementation (standard)
% Options for call to computeMultiPointTrans with block assembly splitted in blocks (necessary for large system)
opts = {'eta'              , eta      , ...
        'verbose'          , true     , ...
        'useTensorAssembly', true};
mpfastruct1 = computeMultiPointTrans(G, rock, opts{:});
% Options for incompMPFA for tensor assembly call
opts = {'useTensorAssembly', true      , ...
        'bc'               , bc        , ...
        'outputFlux'       , true};
state = incompMPFA(state0, G, mpfastruct1, [], opts{:});
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - TA - standard';
caseno         = caseno + 1;

%% mpfa - Tensor Assembly implementation - block assembly (necessary for large systems)
% Options for call to computeMultiPointTrans with block assembly splitted in blocks (necessary for large system)
opts = {'eta'              , eta      , ...
        'blocksize'        , blocksize, ...
        'verbose'          , true     , ...
        'useTensorAssembly', true};
mpfastruct2 = computeMultiPointTrans(G, rock, opts{:});
% Options for incompMPFA for tensor assembly call
opts = {'useTensorAssembly', true      , ...
        'bc'               , bc        , ...
        'outputFlux'       , true};
state = incompMPFA(state0, G, mpfastruct2, [], opts{:});
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - TA - block';
caseno         = caseno + 1;

close all
for i = 1 : numel(vecs)
    figure
    plot(vecs{i}(:, 1), vecs{i}(:, 2), '*-');
    xlabel('z');
    ylabel('pressure');
    title(titles{i});
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
