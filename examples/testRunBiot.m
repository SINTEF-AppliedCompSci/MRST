mrstModule add vemmech mpsaw vem mpfa

clear all
close all

% Switch this on for perturbed grids
pert = 0;
% Test case number (see definitions below)
testCase = 1;
% Constant (between 0 and 1) in MPSA and MPFA which defines position of the continuity point
eta = 1/3;

switch testCase
  case 1
    % Cartesian grid
    gridType = 1;
  case 2
    % Cartesian grid
    gridType = 1;
    eta = 0;
  case 3
    % Triangular grid, 90 degree angles
    gridType = 2;
    eta = 1/3;
  case 4
    % Equilateral triangles
    gridType = 3;
    eta = 1/3;
end

Nx = [30, 30];
G = gridForConvTest(Nx, gridType);

d = G.griddim;

% Physical parameters
mu     = 1;
lambda = 0;
alpha  = 1;
K      = 1;
tau    = 1; % (must be set to 1)
rho    = 1;

params = struct('mu'    , mu    , ...
                'lambda', lambda, ...
                'alpha' , alpha , ...
                'K'     , K     , ...
                'tau'   , tau   , ...
                'rho'   , rho   , ...
                'eta'   , eta);

% Compute analytical solution
output = analyticalBiot(d, params);

params.u_fun      = output.u_fun;
params.p_fun      = output.p_fun;
params.stress_fun = output.stress_fun;
params.force_fun  = output.force_fun;
params.src_fun    = output.src_fun;

% Compute numerical solution
output = runBiotConvSim(G, params);

u = output.u;
p = output.p;

dotest = true;
if dotest
    tbls = output.tbls;
    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);
    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    cellfacetbl = tbls.cellfacetbl;
    bccellfacetbl = crossIndexArray(cellfacetbl, bcfacetbl, {'faces'});
    bccelltbl = projIndexArray(bccellfacetbl, {'cells'});
    bccells = bccelltbl.get('cells');
    bcp = p(bccells);
end    

doplot = true;
if doplot
    
    figure
    plotCellData(G, p);
    title('Pressure, numerical solution');
    
    figure
    plotCellData(G, u(:, 1));
    title('Displacement, x-direction, numerical solution');
    
    figure
    plotCellData(G, u(:, 2));
    title('Displacement, y-direction, numerical solution');

    % prepare input for analytical functions
    for idim = 1 : d
        cc{idim} = G.cells.centroids(:, idim);
    end
    
    figure
    p_exact = params.p_fun(cc{:});
    plotCellData(G, p_exact);
    title('Pressure, analytical solution');
    
    figure
    u1_exact = params.u_fun{1}(cc{:});
    plotCellData(G, u1_exact);
    title('Displacement, x-directiom, analytical solution');

    figure
    u2_exact = params.u_fun{2}(cc{:});
    plotCellData(G, u2_exact);
    title('Displacement, y-directiom, analytical solution');

end


