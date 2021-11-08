function [state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, varargin)
%Set up all inputs to simulateScheduleAD from deck-file
%
% SYNOPSIS:
%   [state0, model, schedule, nonlinear] = initEclipseProblemAD(deck)
%
% REQUIRED PARAMETERS:
%   deck - Representation of the deck-case. Either a string or output from
%          readEclipseDeck from the deckformat module.
%
% OPTIONAL PARAMETERS:
%   useMexGeometry     - Use mcomputeGeometry to compute geometry data.
%                        Default: useMex optional value.
%   useMexProcessGrid  - Use mprocessGRDECL to build grid. Default: useMex
%   TimestepStrategy   - Strategy to pick time-steps in NonLinearSolver
%                        output.
%   useCPR             - Prefer CPR over general iterative solvers for
%                        fully-implicit systems
%   useMex             - Use C/C++ acceleration where possible (e.g. in AD
%                        backend)
%   AutoDiffBackend    - Name of AD backend to use (sparse/diagonal/diagonal-mex)
%   model              - Model to use. Can be used to avoid setup.
%   G                  - Grid. Can be used to avoid setting up grid.
%   getSchedule        - Can be disabled to avoid processing schedule.
%   getInitialState    - Create initial state from deck (EQUIL or direct
%                        assignment of initial conditions)
%   splitDisconnected  - Passed onto processGRDECL.
%   useLegacyModels    - Whether or not to construct the original,
%                        monolithic physical models.  Stop-gap solution
%                        until all examples have been ported to the
%                        Generic* framework and only supported for
%                        three-phase black-oil w/polymer
%                        (`ThreePhaseBlackOilPolymerModel`).
%
% RETURNS:
%   state0    - Initial state.
%   model     - A model instance which attempts to realize all options in
%               the deck.
%   schedule  - The schedule containing the time-steps and driving forces.
%   nonlinear - NonLinearSolver class instance, with attached LinearSolver
%               and timeStepSelector classes suitable for simulating the
%               problem.
% NOTE:
%   name
%
% SEE ALSO:
%   simulateScheduleAD, getNonLinearSolver, packSimulationProblem,
%   selectLinearSolverAD

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

    gravity reset on;

    mrstModule add deckformat ad-core ad-blackoil ad-props
    if ischar(deck)
        % Path to some deck?
        deck = readEclipseDeck(deck);
    end
    deck = convertDeckUnits(deck);

    opt = struct('useMexGeometry',       [], ...
                 'useMexProcessGrid',    [], ...
                 'TimestepStrategy',     'iteration', ...
                 'useCPR',               true, ...
                 'useMex',               mrstSettings('get', 'useMEX'), ...
                 'pvtMethodOil',         'parallel', ...
                 'pvtMethodGas',         'linshift', ...
                 'optimizeTables',       false, ...
                 'rowMajorAD',           false, ...
                 'AutoDiffBackend',      [],    ...
                 'UniformFacilityModel', false, ...
                 'model',                [],    ...
                 'G',                    [],    ...
                 'permTolerance',        0,     ...
                 'maxIterations',        12,    ...
                 'useRelaxation',        true,  ...
                 'deckfn',               [],    ...
                 'getSchedule',          true,  ...
                 'getInitialState',      true,  ...
                 'SplitDisconnected',    false, ...
                 'PreserveCpNodes',      true,  ...
                 'useLegacyModels',      false);
    [opt, extra] = merge_options(opt, varargin{:});
    if ~isempty(opt.deckfn)
        deck = opt.deckfn(deck);
    end
    if isempty(opt.useMexGeometry)
        opt.useMexGeometry = opt.useMex;
    end
    if isempty(opt.useMexProcessGrid)
        opt.useMexProcessGrid = opt.useMex;
    end
    if isempty(opt.model)
        model = initializeModel(deck, opt);
    else
        model = opt.model;
    end
    
    % Set faster backend if grid is sufficiently large and supported by the
    % model.
    if isempty(opt.AutoDiffBackend) 
        if (isa(model, 'NaturalVariablesCompositionalModel') || ...
            isa(model, 'ThreePhaseBlackOilModel')) && model.G.cells.num > 1000
            if opt.useMex
                opt.AutoDiffBackend = 'diagonal-mex';
            else
                opt.AutoDiffBackend = 'diagonal';
            end
        else
            opt.AutoDiffBackend = 'sparse';
        end
    end
    switch lower(opt.AutoDiffBackend)
        case 'diagonal-mex'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', opt.rowMajorAD);
        case 'diagonal'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', opt.useMex, 'rowMajor', opt.rowMajorAD);
        case 'sparse'
            % Do nothing, this is default;
        otherwise
            error('Unknown backend string %s', opt.AutoDiffBackend);
    end
    % If components are not present and all wells are simple, we can use
    % the uniform facility model, which is vectorized and faster when many
    % wells are present
    cnames = model.getComponentNames();
    if isempty(cnames) && opt.UniformFacilityModel && ~isa(model, 'GenericReservoirModel')
        model.FacilityModel = UniformFacilityModel(model);
    end
    % Set up state
    if opt.getInitialState
        state0 = initStateDeck(model, deck);
    else
        state0 = [];
    end
    % Set up schedule
    if opt.getSchedule
        schedule = convertDeckScheduleToMRST(model, deck, extra{:});
    else
        schedule = [];
    end

    if nargout > 3
        nonlinear = getNonLinearSolver(model, ...
                                       'TimestepStrategy', opt.TimestepStrategy, ...
                                       'useCPR',           opt.useCPR,           ...
                                       'maxIterations',    opt.maxIterations,    ...
                                       'useRelaxation',    opt.useRelaxation);
    end
end

function model = initializeModel(deck, opt)
    % Set up grid
    rock  = initEclipseRock(deck);
    
    if isempty(opt.G)
        if ~isfield(deck.GRID, 'ACTNUM')
            nc = prod(deck.GRID.cartDims);
            deck.GRID.ACTNUM = ones(nc, 1);
        end
            
        if isfield(rock, 'ntg')
            pv = rock.poro.*rock.ntg;
        else
            pv = rock.poro;
        end
        perm_ok = any(rock.perm > opt.permTolerance, 2);
        deck.GRID.ACTNUM = double(deck.GRID.ACTNUM > 0 & pv > 0 & perm_ok);

        G = initEclipseGrid(deck, 'SplitDisconnected', opt.SplitDisconnected, ...
                                  'useMex', opt.useMexProcessGrid, ...
                                  'PreserveCpNodes', ...
                                  ~opt.useMexProcessGrid && opt.PreserveCpNodes);
        if numel(G) > 1
            warning('Multiple disconnected grids found. Picking largest.');
            G = G(1);
        end
    else
        G = opt.G;
    end
    if ~isfield(G.cells, 'centroids')
        if opt.useMexGeometry
            mrstModule add libgeometry
            G = mcomputeGeometry(G);
        else
            G = computeGeometry(G);
        end
    end
    fluid = initDeckADIFluid(deck, 'G', G, 'useMex', opt.useMex, ...
                                   'optimize',     opt.optimizeTables, ...
                                   'pvtMethodGas', opt.pvtMethodGas, ...
                                   'pvtMethodOil', opt.pvtMethodOil);
    
    rock  = compressRock(rock, G.cells.indexMap);
    model = selectModelFromDeck(G, rock, fluid, deck, ...
                                'UseLegacyModels', opt.useLegacyModels);
    model.dpMaxRel = 0.2;
end
