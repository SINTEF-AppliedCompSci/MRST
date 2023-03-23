function [G_global,TPFAoperators,fracplanes] = EDFMpreprocessor(G_matrix, fracplanes, varargin)
% preProcessingFractures identifies fracture-matrix connections, imposes a
% fracture grid and computes the fracture-matrix conductivity index.
%
% SYNOPSIS:
%   [G,fracplanes] = preProcessingFractures(G, fracplanes)
%   [G,fracplanes] = preProcessingFractures(G, fracplanes, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G_matrix          - Grid data structure containing geometrical
%                       information for matrix.
%
%   fracplanes  - 1-by-n structure where n is the number of fractures. Each
%                 column contains a set of coplanar points that define
%                 the fracture polygon, a value for the fracture
%                 aperture and fracture porosity and permeability values.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
%
%   fractureCellSize   - Dimensionless element size (>0 and <1) for the
%                        fracture grid.
%
%   Tolerance    - numerical tolerance. 
%                        Default: 1e-5
%
% RETURNS:
%   G_global - Global grid structure containing both matrix and fracture cells and
%       information about fracture-matrix NNC connections
%
%   TPFAoperators - operators to replace model.operators with. After
%                   generating reservoir model using
%                   ThreePhaseBlackOilModel, use the command line
%                   model.operators=TPFAoperators.
%
%   fracplanes  - Structure with added
%
% SEE ALSO:
%   markcells, getPlaneNormals, gridPlanarFracture, fracMatrixConnections


opt = struct('fractureCellSize'  ,    -1     , ...
             'Tolerance'   ,    1e-5, ...
             'plotgrid', false);
         
opt = merge_options(opt, varargin{:});

tol=opt.Tolerance;

fracplanes = getPlaneNormals(fracplanes);


%% Generate fracture grids (use gridPlanarFracture_mod)
% Create FracGrid which contains fracture grids for every fracture in
% fracplanes. FracGrid.Frac1,....,N will have a grid structure just like
% the matrix grid. Additionally, they will also contain global starting
% indices for the cells, nodes and faces, e.g. if starting global index for
% Frac1 is 1001, then the second cell in Frac1 will be the 1002-th cell in
% the global grid. The porosity and permeability data from fracplanes will
% be transferred into FracGrid. 
FracGrid=struct;
cstart = G_matrix.cells.num + 1;
fstart = G_matrix.faces.num + 1;
nstart = G_matrix.nodes.num + 1;

for i = 1:numel(fracplanes)
    fieldname = ['Frac',num2str(i)];
%     isrectangular = ismember(i,rectangular);
    K_frac=fracplanes(i).perm;
    phi_frac=fracplanes(i).poro;
    FracGrid.(fieldname) = gridPlanarFracture_mod(G_matrix, fracplanes(i), tol,'cellSize',opt.fractureCellSize);
    FracGrid.(fieldname).cells.start = cstart;
    FracGrid.(fieldname).faces.start = fstart;
    FracGrid.(fieldname).nodes.start = nstart;
    FracGrid.(fieldname).rock.perm = ones(FracGrid.(fieldname).cells.num,1)*K_frac;
    FracGrid.(fieldname).rock.poro = ones(FracGrid.(fieldname).cells.num,1)*phi_frac;
    
    cstart = cstart + FracGrid.(fieldname).cells.num;
    fstart = fstart + FracGrid.(fieldname).faces.num;
    nstart = nstart + FracGrid.(fieldname).nodes.num;   
end

%% Plot Grid
if opt.plotgrid
    figure; plotGrid(G_matrix,'facealpha',0);
    for i = 1:numel(fieldnames(FracGrid))
        plotGrid(FracGrid.(['Frac',num2str(i)]));
    end
    view(15,20);
end

%% Generate global grid. Fracture and matrix grids are saved under
% G_global.FracGrid and G_global.Matrix (use assembleGlobalGrid)
G_global=G_matrix; G_global.FracGrid=FracGrid; G_global.nnc=[]; % this line is a trick to fit the input for assembleGlobalGrid
G_global=assembleGlobalGrid(G_global);
% at this point, G_global does not contain any nnc information. It will
% contain G_global.Matrix and G_global.FracGrid which will each contain the
% original matrix and fracture grids.

%% Perform intersection check between matrix and fracture grids. 
% Assign nnc's, areas, CI's to G_global.nnc.cells, G_global.nnc.CI and G_global.nnc.area. For
% G_global.nnc.cells, the indices assigned are the global grid cell indices.
%   -assume all intersection areas are convex, use convexpolygonarea to
%   calculate areas
%   -use computeEffectiveTrans to generate G_global.nnc.T from G_global.nnc.CI
G_global=fracturematrixNNC3D(G_global,tol);
G_global=computeEffectiveTrans(G_global);

%% Perform intersection check between fracture and fracture grids. Assign
% connectivity and transmissibilities to G_global.nnc.cells and G_global.nnc.T. Again,
% the cell indices are global indices. Fracplanes will have new field that
% indicates which fractures intersect with it.
% (use fracturefractureNNCs3D)
[G_global,fracplanes]=fracturefractureNNCs3D(G_global,fracplanes,tol);

%% Run the following to generate neighbouring transmissibilities
T = computeTrans(G_global, G_global.rock);
cf = G_global.cells.faces(:,1);
nf = G_global.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G_global.nnc.T];

%% Run the following to generate TPFA operators for ThreePhaseBlackOilModel
N = getNeighbourship(G_global, 'topological', true);
intx = all(N ~= 0, 2);

% Send in internal transmissibility and neighborship to the operator setup
TPFAoperators = setupOperatorsTPFA(G_global, G_global.rock, 'trans', T(intx), 'neighbors', N(intx, :));

end