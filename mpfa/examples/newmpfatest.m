
mrstModule add mimetic mpfa incomp

isverbose = false;
eta       = 1/3;
blocksize = 1000;

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 100; ny = 100;
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

titles = {'mpfa - jostein', 'mpfa - standard', 'mpfa - block', 'mpfa - well', 'mpfa - well - block'};

clear mpfastructs states pressures fluxes
caseno = 1;

% mpfa - jostein
tic
T_mpfa = computeMultiPointTrans(G, rock, 'eta', eta);
texec = toc;
states{caseno} = initResSol(G, 0, 1);
states{caseno} = incompMPFA(states{caseno}, G, T_mpfa, fluid, 'wells', W);
pressures{caseno} = states{caseno}.pressure;
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
titles{caseno} = 'mpfa - jostein';
caseno = caseno + 1;

% mpfa - standard
tic
mpfastructs{caseno} = computeMultiPointTrans2(G, rock, 'eta', eta, 'verbose', ...
                                              isverbose);
texec = toc;
states{caseno} = incompMPFA2(G, mpfastructs{caseno}, 'wells', W);
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno} = states{caseno}.flux;
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
titles{caseno} = 'mpfa - standard';
caseno = caseno + 1;


% mpfa - block
tic
mpfastructs{caseno} = computeMultiPointTrans2(G, rock, 'eta', eta, 'blocksize', ...
                                              blocksize, 'verbose', isverbose);
texec = toc;
states{caseno} = incompMPFA2(G, mpfastructs{caseno}, 'wells', W);
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno} = states{caseno}.flux;
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
titles{caseno} = 'mpfa - block';
caseno = caseno + 1;

% mpfa - well
tic
mpfastructs{caseno} = computeNeumannMultiPointTrans(G, rock, 'eta', eta, 'verbose', isverbose);
texec = toc;
states{caseno} = incompMPFA3(G, mpfastructs{caseno}, W, 'outputFlux', true);
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno} = states{caseno}.flux;
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
titles{caseno} = 'mpfa - well';
caseno = caseno + 1;

% mpfa - well - block
tic
mpfastructs{caseno} = computeNeumannMultiPointTrans(G, rock, 'eta', 1/3, 'verbose', isverbose);
texec = toc;
states{caseno} = incompMPFA3(G, mpfastructs{caseno}, W, 'outputFlux', true);
pressures{caseno} = states{caseno}.pressure;
fluxes{caseno} = states{caseno}.flux;
fprintf('Done with %s in %g sec\n', titles{caseno}, texec);
titles{caseno} = 'mpfa - well - block';
caseno = caseno + 1;



%% Plotting
for i = 1 : numel(pressures)
    figure(i)
    clf
    plotCellData(G, pressures{i} ./ barsa());
    title(titles{i});
    colorbar
end

return

%% check mass conservation
extfaces = any(G.faces.neighbors == 0, 2);
intfacetbl.faces = find(~extfaces);
intfacetbl.num   = numel(intfacetbl.faces);
tbls = mpfastruct.tbls;
cellfacetbl = tbls.cellfacetbl;
cfc = cellfacetbl.cells;
cff = cellfacetbl.faces;
sgn = 2*(cfc == G.faces.neighbors(cff, 1)) - 1;
op1 = setupTableMapping(intfacetbl, cellfacetbl, {'faces'});
op1 = sparse(1 : cellfacetbl.num, 1 : cellfacetbl.num, sgn)*op1;
celltbl.cells = (1 : G.cells.num)';
celltbl.num   = G.cells.num;
op2 = setupTableMapping(cellfacetbl, celltbl, {'cells'});
div = op2*op1;

massbal = div*flux;

%% plotting

figure(2)
clf
plotCellData(G, massbal);
title('Mass balance'); view(2), axis tight off
colorbar('Location','SouthOutside');

%% solve for mimetic

IP = computeMimeticIP(G, rock, 'InnerProduct', 'ip_simple');
resSol2 = initState(G, W, 0);

%% Solve mimetic linear hybrid system
resSol2 = incompMimetic(resSol2, G, IP, fluid, 'wells', W);
pressure = resSol2.pressure;

figure(3)
clf
plotCellData(G, pressure ./ barsa());
title('Pressure: mimetic'); view(2), axis tight off
colorbar('Location','SouthOutside');