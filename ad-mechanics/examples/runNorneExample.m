function [model, states, initState] = runNorneExample(varargin)
%% Example: Poroelasticity simulation applied to the Norne case.
% 
% The simulation options are gathered in the opt structure. If opt=[] the
% simulation is run with the default options defined below
%
% **  Summary of the options ** 
%
% option 'norne_case' :
%
%     * 'full'       : 7392 cells
%     * 'mini Norne' :  605 cells
%
% option 'bc_case' :
%
%     * 'no displacement' : All nodes belonging to external faces have displacement
%                           equal to zero
%     * 'bottom fixed'    : The nodes that belong to the bottom have zero
%                           displacement, while a given pressure is imposed on
%                           the external faces that are not bottom faces.
%
% option 'method' :
%
%     * 'fully coupled'          : The mechanical and flow equations are solved fully coupled.
%     * 'fixed stress splitting' : The mechanical and flow equations are solved
%                                  sequentially using a fixed stress splitting
%
% option 'fluid_model' :
%
%     * 'blackoil'  : blackoil model is used for the fluid (gas is injected, see
%                     schedule below)
%     * 'oil water' : Two phase oil-water
%     * 'water'     : water model is used for the fluid

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

    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

    % setup default option values
    opt = struct('norne_case'         , 'full'          , ...
                 'bc_case'            , 'bottom fixed'  , ...
                 'method'             , 'fully coupled' , ...
                 'fluid_model'        , 'water'         , ...
                 'nonlinearTolerance' , 1e-6            , ...
                 'splittingTolerance' , 1e-3            , ...
                 'verbose'            , false           , ...
                 'splittingVerbose'   , false);
    opt = merge_options(opt, varargin{:});

    % overwrite the default options by the given option
    % optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
    % optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
    % opt = merge_options(default_opt, optlist{:});


    %% Load Norne grid

    if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL()),
        error('Unable to obtain simulation model subset');
    end


    grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
    grdecl = readGRDECL(grdecl);
    fullCartDims = grdecl.cartDims;
    usys   = getUnitSystem('METRIC');
    grdecl = convertInputUnits(grdecl, usys);
    switch opt.norne_case
      case 'full'
        grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]); 
      case 'full padded'
        grdeclo = cutGrdecl(grdecl, [10 25; 35 55; 1 22]); 
        oncells = prod(grdeclo.cartDims);
        oindex = 1:oncells; 
        oindex = reshape(oindex, grdeclo.cartDims); 
        grdecl = padGrdecl(grdeclo, [true, true, true], [[60 50; 40 40] * 10; 10 10], 'relative', true); 
        grdecl.ACTNUM = ones(prod(grdecl.cartDims), 1); 
        indexmapn = zeros(grdecl.cartDims); 
        indexmapn(2:end - 1, 2:end - 1, 2:end - 1) = oindex; 
        indexmapn(indexmapn == 0) = prod(grdeclo.cartDims) + 1; 
        ocells = indexmapn < prod(grdeclo.cartDims);
        oindex = indexmapn;
      case 'mini Norne'
        grdecl = cutGrdecl(grdecl, [10 20; 35 45; 1 5]); 
      otherwise
        error('norne_case not recognized');
    end
    grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
    G = processGRDECL(grdecl);
    G = G(1);
    G = computeGeometry(G);

    %% Setup rock parameters (for flow)

    switch opt.norne_case
      case {'full', 'mini Norne'}
        perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
        rock.perm = perm(G.cells.indexMap, :);
        rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);
      case 'full padded'
        perm = [grdeclo.PERMX, grdeclo.PERMY, grdeclo.PERMZ;
                [1 1 1]*milli*darcy ];
        poro = [grdeclo.PORO;0.1];    
        rock.perm = perm(oindex(G.cells.indexMap), :);
        rock.poro = max(poro(oindex(G.cells.indexMap)), 0.1);
        [i, j, k] = ind2sub(G.cartDims, G.cells.indexMap); 
        ind = (k>1 & k<G.cartDims(3) &  oindex(G.cells.indexMap) == numel(poro)); 
        rock.perm(ind, :) = rock.perm(ind, :) * 100;
      otherwise
        error('norne_case not recognized');
    end


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
                                   1000*kilogram/meter^3, 'c', 1e-10, ...
                                   'cR', 4e-10, 'pRef', pRef);
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

        ind = (G.faces.neighbors(:, 1) == 0 | G.faces.neighbors(:, 2) == 0);
        ind = find(ind);
        nodesind = mcolon(G.faces.nodePos(ind), G.faces.nodePos(ind + 1) - 1);
        nodes = G.faces.nodes(nodesind);
        bcnodes = zeros(G.nodes.num);
        bcnodes(nodes) = 1;
        bcnodes = find(bcnodes == 1);
        nn = numel(bcnodes);
        u = zeros(nn, 3);
        m = ones(nn, 3);
        disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);
        force_bc = [];

      case 'bottom fixed'

        nx = G.cartDims(1);
        ny = G.cartDims(2);
        nz = G.cartDims(3);

        % Find the bottom nodes. On these nodes, we impose zero displacement

        c = zeros(nx*ny*nz, 1);
        c(G.cells.indexMap) = (1 : numel(G.cells.indexMap))';
        bottomcells = c(nx*ny*(nz - 1) +  (1 : (nx*ny))');
        faces = G.cells.faces(mcolon(G.cells.facePos(bottomcells), G.cells.facePos(bottomcells ...
                                                          + 1) - 1), :);
        bottomfaces = faces( faces(:, 2) == 6  , 1);
        indbottom_nodes = mcolon(G.faces.nodePos(bottomfaces), G.faces.nodePos(bottomfaces ...
                                                          + 1) - 1);
        bottom_nodes = G.faces.nodes(indbottom_nodes);
        isbottom_node = false(G.nodes.num, 1);
        isbottom_node(bottom_nodes) = true;
        bcnodes = find(isbottom_node);

        nn = numel(bcnodes);
        u = zeros(nn, 3);
        m = ones(nn, 3);
        disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);

        % Find outer faces that are not at the bottom. On these faces, we impose
        % a given pressure.

        is_outerface1 = (G.faces.neighbors(:, 1) == 0);
        is_outerface1(bottomfaces) = false;
        is_outerface2 = G.faces.neighbors(:, 2) == 0;
        is_outerface2(bottomfaces) = false;

        is_outerface = is_outerface1 | is_outerface2;

        outer_faces = find(is_outerface);

        outer_pressure = pRef;
        signcoef = (G.faces.neighbors(outer_faces, 1) == 0) - ...
            (G.faces.neighbors(outer_faces, 2) == 0);
        n = bsxfun(@times, G.faces.normals(outer_faces, :), signcoef./ ...
                   G.faces.areas(outer_faces));
        force = bsxfun(@times, n, outer_pressure);

        force_bc = struct('faces', outer_faces, 'force', force);


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
    gravity on;


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
    refdepth = G.cells.centroids(1, 3); % for example...
    injcell  = 10; % for example...
    prodcell = G.cells.num;
    
    if strcmp(opt.norne_case, 'full padded')
        injcello = 10; % for example
        pos = [0.457517456153883   7.321672863418036   0.002564025720821] * 1e6; 
        [~, injcell] = min(sum(bsxfun(@minus, G.cells.centroids, pos).^2, 2))
        injcell = find(oindex(:) == injcello);
        prodcello = 7392; % for example
        pos = [0.458944373173747   7.322630019653659   0.002783994408426] * 1e6; 
        [~, prodcell] = min(sum(bsxfun(@minus, G.cells.centroids, pos).^2, 2)); 
        prodcell = find(oindex(:) == prodcello);
    end
    
    W = addWell(W, G, rock, injcell, ...
                'Type'    , 'rate', ...
                'Val'     , 2.5e6/day, ...
                'Sign'    , 1,  ...
                'Comp_i'   , [0, 0, 1], ... % inject gas
                'Name'    , 'inj',  ...
                'refDepth', refdepth);

    
    W = addWell(W, G, rock, prodcell, ...
                'Type'    ,'bhp', ...
                'Val'     , pRef, ...
                'Sign'    , -1,  ...
                'Comp_i'   , [0, 1, 0], ... % one-phase test case
                'Name'    , 'prod',  ...
                'refDepth', refdepth);

    switch opt.fluid_model
      case 'blackoil'
        W(1).compi = [0, 0, 1];
        W(2).compi = [0, 1, 0];
      case 'oil water'
        W(1).compi = [1 0];
        W(1).val   = 1e4/day;
        W(2).compi = [0 1];
      case 'water'
        W(1).compi = [1];
        W(1).val   = 1e4/day;
        W(2).compi = [1];
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
    initState   = computeInitDisp(model, initState, [], 'pressure', initState.pressure);

    [wellSols, states, schedulereport] = simulateScheduleAD(initState, model, schedule);

end
