function [state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, varargin)
    mrstModule add deckformat ad-core ad-blackoil ad-props
    deck = convertDeckUnits(deck);

    opt = struct('useMexGeometry', false, ...
                 'TimestepStrategy', 'iteration', ...
                 'useMex',          false, ...
                 'AutoDiffBackend',  [], ...
                 'UniformFacilityModel', false, ...
                 'model',           [], ...
                 'G',               [], ...
                 'getSchedule',     true, ...
                 'getInitialState', true, ...
                 'SplitDisconnected', false);
    [opt, extra] = merge_options(opt, varargin{:});
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
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
        case 'diagonal'
            model.AutoDiffBackend = DiagonalAutoDiffBackend();
        case 'sparse'
            % Do nothing, this is default;
    end
    % If components are not present and all wells are simple, we can use
    % the uniform facility model, which is vectorized and faster when many
    % wells are present
    cnames = model.getComponentNames();
    if isempty(cnames) && opt.UniformFacilityModel && ~isa(model, 'ExtendedReservoirModel')
        model.FacilityModel = UniformFacilityModel(model);
    end
    % Set up schedule
    if opt.getSchedule
        schedule = convertDeckScheduleToMRST(model, deck, extra{:});
    else
        schedule = [];
    end
    % Set up state
    if opt.getInitialState
        state0 = initStateDeck(model, deck);
    else
        state0 = [];
    end
    
    if nargout > 3
        nonlinear = NonLinearSolver();
        nonlinear.LinearSolver = selectLinearSolverAD(model);
        switch lower(opt.TimestepStrategy)
            case 'none'
                % Do nothing
                sel = [];
            case 'iteration'
                % Control on iterations
                sel = IterationCountTimeStepSelector('targetIterationCount', 8);
            case 'ds'
                % Control on saturation change
                sel = ...
                    StateChangeTimeStepSelector('targetProps', {'s'},...
                                                'targetChangeAbs', 0.2, ...
                                                'targetIterationCount', inf);
            case 'dsdc'
                % Control on saturation + components
                names = {'s'};
                targets = 0.2;
                if isa(model, 'ThreePhaseCompositionalModel')
                    names = {'s', 'components'};
                    targets = [0.2, 0.2];
                end
                sel = ...
                    StateChangeTimeStepSelector('targetProps', names,...
                                                'targetChangeAbs', targets, ...
                                                'targetIterationCount', inf);

            otherwise
                error('Unknown timestepping strategy %s', opt.TimestepStrategy);
        end
        if ~isempty(sel)
            sel.firstRampupStepRelative = 0.1;
            sel.firstRampupStep = 1*day;
            nonlinear.timeStepSelector = sel;
        end
    end
end

function model = initializeModel(deck, opt)
    % Set up grid
    rock  = initEclipseRock(deck);
    
    if isempty(opt.G)
        if isfield(deck.GRID, 'ACTNUM')
            if isfield(rock, 'ntg')
                pv = rock.poro.*rock.ntg;
            else
                pv = rock.poro;
            end
            perm_ok = ~all(rock.perm == 0, 2);
            deck.GRID.ACTNUM = double(deck.GRID.ACTNUM > 0 & pv > 0 & perm_ok);
        end

        G = initEclipseGrid(deck, 'SplitDisconnected', opt.SplitDisconnected);
        if numel(G) > 1
            warning('Multiple disconnected grids found. Picking largest.');
            G = G(1);
        end
    else
        G = opt.G;
    end
    if opt.useMexGeometry
        mrstModule add libgeometry
        G = mcomputeGeometry(G);
    else
        G = computeGeometry(G);
    end
    fluid = initDeckADIFluid(deck, 'G', G);
    gravity reset on;
    
    rock  = compressRock(rock, G.cells.indexMap);
    model = selectModelFromDeck(G, rock, fluid, deck);
    model.dpMaxRel = 0.2;
end
