function [model, states] = run2DCase(varargin)
%
%
% SYNOPSIS:
%   function run2DCase(varargin)
%
% DESCRIPTION:
%    Example which setups poroelasticity computations for a two dimensional domain.
%
%    Simple well setup: One injection well in the middle.
%
%
% PARAMETERS:
%   varargin - see options below
%
% SEE ALSO: runAllNorneExamples

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

%% 2D example for poroelasticity
% Options that can be chosen in this example (see opt structure below)
%
% - fluid model :
%
%   * 'water'     : single phase model
%   * 'oil water' : two phases model
%   * 'blackoil'  : Three phases, blackoil model
%
% - solver :
%
%   * 'fully coupled'          : fully coupled solver
%   * 'fixed stress splitting' : solver using fixed stress splitting
%
%
% - Cartesian grid :
%
%   * 'cartDim' : number of cells in x and y directions
%   * 'L'       : Physical length in x and y directions
%
% - boundary conditions for the mechanics (only one choice here)
%
%   * 'bottom fixed' : Fixed bottom



    %% Load required modules

    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

    %% Setup default options
    opt = struct('cartDim'            , [100, 10], ...
                 'L'                  , [100, 10], ...
                 'fluid_model'        , 'water', ...
                 'method'             , 'fully coupled', ...
                 'bc_case'            , 'bottom fixed', ...
                 'nonlinearTolerance' , 1e-6, ...
                 'splittingTolerance' , 1e-6, ...
                 'verbose'            , false, ...
                 'splittingVerbose'   , false);

    opt = merge_options(opt, varargin{:});
    %% Setup grid


    G = cartGrid(opt.cartDim, opt.L);
    G = computeGeometry(G);

    %% Setup rock parameters (for flow)

    rock.perm = darcy*ones(G.cells.num, 1);
    rock.poro = 0.3*ones(G.cells.num, 1);


    %% Setup fluid parameters from SPE1

    pRef = 270*barsa;
    switch opt.fluid_model
      case 'blackoil'
        pth = getDatasetPath('spe1');
        fn  = fullfile(pth, 'BENCH_SPE1.DATA');
        deck = readEclipseDeck(fn);
        deck = convertDeckUnits(deck);
        fluid = initDeckADIFluid(deck);
        if isfield(fluid, 'pcOW')
            fluid = rmfield(fluid, 'pcOW');
        end
        if isfield(fluid, 'pcOG')
            fluid = rmfield(fluid, 'pcOG');
        end

        % Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
        fluid.krW = @(s) s.^2;
        fluid.krG = @(s) s.^2;
        fluid.krOW = @(s) s.^2;
        fluid.krOG = @(s) s.^2;
        pRef = deck.PROPS.PVTW(1);

      case {'oil water'}
        fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 10]*centi*poise, ...
                                   'n',  [1, 1], 'rho', [1000, 700]*kilogram/ ...
                                   meter^3, 'c', 1e-10*[1, 1], 'cR', 4e-10, ...
                                   'pRef', pRef);

      case {'water'}
        fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                                   1000*kilogram/meter^3, 'c', 1e-10, 'cR', ...
                                   4e-10, 'pRef', pRef);
      otherwise
        error('fluid_model  not recognized.');
    end


    %% Setup material parameters for Biot and mechanics

    E          = 1 * giga * Pascal; % Young's module
    nu         = 0.3;               % Poisson's ratio
    alpha      = 1;                 % Biot's coefficient

    % Transform these global properties (uniform) to cell values.
    E          = repmat(E, G.cells.num, 1);
    nu         = repmat(nu, G.cells.num, 1);
    rock.alpha = repmat(alpha, G.cells.num, 1);


    %% Setup boundary conditions for mechanics (no displacement)

    switch opt.bc_case

      case 'no displacement'
        error('not implemented yet');
        ind = (G.faces.neighbors(:, 1) == 0 | G.faces.neighbors(:, 2) == 0);
        ind = find(ind);
        nodesind = mcolon(G.faces.nodePos(ind), G.faces.nodePos(ind + 1) - 1);
        nodes = G.faces.nodes(nodesind);
        bcnodes = zeros(G.nodes.num);
        bcnodes(nodes) = 1;
        bcnodes = find(bcnodes == 1);
        nn = numel(bcnodes);
        u = zeros(nn, 2);
        m = ones(nn, 2);
        disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);
        force_bc = [];

      case 'bottom fixed'

        nx = G.cartDims(1);
        ny = G.cartDims(2);

        % Find the bottom nodes. On these nodes, we impose zero displacement

        c = zeros(prod(G.cartDims), 1);
        c(G.cells.indexMap) = (1 : numel(G.cells.indexMap))';

        bc = pside([], G, 'Ymin', 100);
        bottomfaces = bc.face;
        indbottom_nodes = mcolon(G.faces.nodePos(bottomfaces), ...
                                 G.faces.nodePos(bottomfaces + 1) - 1);
        bottom_nodes = G.faces.nodes(indbottom_nodes);
        bottom_nodes = unique(bottom_nodes);

        nn = numel(bottom_nodes);
        u = zeros(nn, G.griddim);
        m = ones(nn, G.griddim);
        disp_bc = struct('nodes', bottom_nodes, 'uu', u, 'mask', m);

        % Find outer faces that are not at the bottom. On these faces, we impose
        % a given pressure.

        bc = pside([], G, 'Xmin', 100);
        bc = pside(bc, G, 'Xmax', 100);
        bc = pside(bc, G, 'Ymax', 100);
        sidefaces = bc.face;

        signcoef = (G.faces.neighbors(sidefaces, 1) == 0) - (G.faces.neighbors(sidefaces, ...
                                                          2) == 0);
        n = bsxfun(@times, G.faces.normals(sidefaces, :), signcoef./ ...
                   G.faces.areas(sidefaces));
        force = bsxfun(@times, n, pRef);

        force_bc = struct('faces', sidefaces, 'force', force);


      otherwise
        error('bc_cases not recognized')
    end

    el_bc = struct('disp_bc' , disp_bc, ...
                   'force_bc', force_bc);


    %% Setup load for mechanics

    % In this example we do not impose any volumetric force
    loadfun = @(x) (0*x);


    %% Gather all the mechanical parameters in a struct

    mech = struct('E', E, 'nu', nu, 'el_bc', el_bc, 'load', loadfun);


    %% Gravity
    % The gravity in this option affects only the fluid behavior
    gravity off;


    %% Setup model

    modeltype = [opt.method, ' and ', opt.fluid_model];
    fullycoupledOptions = {'verbose', opt.verbose};
    splittingOptions = {'splittingTolerance', opt.splittingTolerance, ...
                        'splittingVerbose', opt.splittingVerbose};
    switch modeltype

      case 'fully coupled and blackoil'
        model = MechBlackOilModel(G, rock, fluid, mech, fullycoupledOptions{: ...
                   });

      case 'fixed stress splitting and blackoil'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'blackoil', ...
                                               splittingOptions{:});

      case 'fully coupled and oil water'
        model = MechOilWaterModel(G, rock, fluid, mech, fullycoupledOptions{: ...
                   });

      case 'fixed stress splitting and oil water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'oil water', ...
                                               splittingOptions{:});

      case 'fully coupled and water'
        model = MechWaterModel(G, rock, fluid, mech, fullycoupledOptions{: });

      case 'fixed stress splitting and water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'water', ...
                                               splittingOptions{:});

      otherwise
        error('modeltype not recognized.');
    end



    %% Setup wells
    W = [];
    refdepth = min(G.cells.centroids(:, G.griddim)) ;
    
    ind = ceil(G.cartDims/2);
    injcell = sub2ind(G.cartDims, ind(1), ind(2));

    warning('off', 'RefDepth:BelowTopConnection');
    
    W = addWell(W, G, rock, injcell, ...
                'Type'    , 'rate', ...
                'Val'     , 1e2/day, ...
                'Sign'    , 1,  ...
                'Comp_i'  , [0, 0, 1], ... % inject gas
                'Name'    , 'inj',  ...
                'refDepth', refdepth);

    % W = addWell(W, G, rock, prodcell, ...
    %             'Type'    ,'bhp', ...
    %             'Val'     , pRef, ...
    %             'Sign'    , -1,  ...
    %             'Comp_i'  , [0, 1, 0], ... % one-phase test case
    %             'Name'    , 'prod',  ...
    %             'refDepth', refdepth);

    switch opt.fluid_model
      case 'blackoil'
        W(1).compi = [0, 0, 1];
      case 'oil water'
        W(1).compi = [1 0];
        W(1).val   = 1e4/day;
      case 'water'
        W(1).compi = 1;
        W(1).val  = 1e-3/day;
      otherwise
        error('fluid_model not recognized.')
    end

    facilityModel = FacilityModel(model.fluidModel);
    facilityModel = facilityModel.setupWells(W);
    model.FacilityModel = facilityModel;



    %% Setup schedule
    schedule.step.val     = [1*day*ones(1, 1); 5*day*ones(20, 1)];
    schedule.step.control = ones(numel(schedule.step.val), 1);
    schedule.control      = struct('W', W);

    %% Setup initial state
    clear initState;
    initState.pressure = pRef*ones(G.cells.num, 1);
    switch opt.fluid_model
      case 'blackoil'
        init_sat = [0, 1, 0];
        initState.rs  = 0.5*fluid.rsSat(initState.pressure);
      case 'oil water'
        init_sat = [0, 1];
      case 'water'
        init_sat = [1];
      otherwise
        error('fluid_model not recognized.')
    end
    initState.s = ones(G.cells.num, 1)*init_sat;
    initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
    initState = addDerivedQuantities(model.mechModel, initState);

    solver = NonLinearSolver('maxIterations', 100);
    [wsol, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', ...
                                        solver);


end
