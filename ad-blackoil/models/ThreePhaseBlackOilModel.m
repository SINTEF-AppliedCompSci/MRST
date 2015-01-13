classdef ThreePhaseBlackOilModel < ReservoirModel
% Three phase with optional dissolved gas and vaporized oil
properties
    % Flag deciding if gas can be dissolved into the oil phase
    disgas
    % Flag deciding if oil can be vaporized into the gas phase
    vapoil

    % Maximum relative Rs/Rv increment
    drsMaxRel
    % Maximum absolute Rs/Rv increment
    drsMaxAbs
end

methods
    function model = ThreePhaseBlackOilModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});

        % Typical black oil is disgas / dead oil, but all combinations
        % are supported
        model.vapoil = false;
        model.disgas = false;

        % Max increments
        model.drsMaxAbs = inf;
        model.drsMaxRel = inf;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.oil = true;
        model.gas = true;
        model.water = true;
        model.saturationVarNames = {'sw', 'so', 'sg'};

        model = merge_options(model, varargin{:});

        d = model.inputdata;
        if ~isempty(d)
            % Assume ECL-style input deck, as this is the only
            % supported format at the moment.
            if isfield(d, 'RUNSPEC')
                if isfield(d.RUNSPEC, 'VAPOIL')
                    model.vapoil = d.RUNSPEC.VAPOIL;
                end
                if isfield(d.RUNSPEC, 'DISGAS')
                    model.disgas = d.RUNSPEC.DISGAS;
                end
            else
                error('Unknown dataset format!')
            end
        end
    end
    
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        switch(lower(name))
            case {'rs', 'rv'}
                fn = lower(name);
                index = 1;
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsBlackOil(state0, state, model, dt, ...
                        drivingForces, varargin{:});

    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        saturations = lower(model.saturationVarNames);
        wi = strcmpi(saturations, 'sw');
        oi = strcmpi(saturations, 'so');
        gi = strcmpi(saturations, 'sg');

        vars = problem.primaryVariables;
        removed = false(size(vars));
        if model.disgas || model.vapoil
            % The VO model is a bit complicated, handle this part
            % explicitly.
            state0 = state;

            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);

            [vars, ix] = model.stripVars(vars, 'pressure');
            removed(~removed) = removed(~removed) | ix;

            % Black oil with dissolution
            so = model.getProp(state, 'so');
            sw = model.getProp(state, 'sw');
            sg = model.getProp(state, 'sg');

            % Magic status flag, see inside for doc
            st = getCellStatusVO(model, state0, so, sw, sg);

            dr = model.getIncrement(dx, problem, 'x');
            dsw = model.getIncrement(dx, problem, 'sw');
            % Interpretation of "gas" phase varies from cell to cell, remove
            % everything that isn't sG updates
            dsg = st{3}.*dr - st{2}.*dsw;

            if model.disgas
                state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                       'rs', model.drsMaxRel, model.drsMaxAbs);
            end

            if model.vapoil
                state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                       'rv', model.drsMaxRel, model.drsMaxAbs);
            end

            dso = -(dsg + dsw);

            ds = zeros(numel(so), numel(saturations));
            ds(:, wi) = dsw;
            ds(:, oi) = dso;
            ds(:, gi) = dsg;

            state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);
            % We should *NOT* be solving for oil saturation for this to make sense
            assert(~any(strcmpi(vars, 'so')));
            state = computeFlashBlackOil(state, state0, model, st);
            state.s  = bsxfun(@rdivide, state.s, sum(state.s, 2));

            %  We have explicitly dealt with rs/rv properties, remove from list
            %  meant for autoupdate.
            [vars, ix] = model.stripVars(vars, {'sw', 'so', 'sg', 'rs', 'rv', 'x'});
            removed(~removed) = removed(~removed) | ix;

        end

        % We may have solved for a bunch of variables already if we had
        % disgas / vapoil enabled, so we remove these from the
        % increment and the linearized problem before passing them onto
        % the generic reservoir update function.
        problem.primaryVariables = vars;
        dx(removed) = [];

        % Parent class handles almost everything for us
        [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);

        % Handle the directly assigned values (i.e. can be deduced directly from
        % the well controls. This is black oil specific.
        W = drivingForces.Wells;
        state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);
    end
end
end
