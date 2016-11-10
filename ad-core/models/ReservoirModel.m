classdef ReservoirModel < PhysicalModel
%Base class for physical models
%
% SYNOPSIS:
%   model = ReservoirModel(G, rock, fluid)
%
% DESCRIPTION:
%   Extension of PhysicalModel class to accomodate reservoir-specific
%   features such as fluid and rock as well as commonly used phases and
%   variables.
%
% REQUIRED PARAMETERS:
%   G     - Simulation grid.
%
%   rock  - Valid rock used for the model.
%
%   fluid - Fluid model used for the model.
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, PhysicalModel

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

properties
    % The fluid model
    fluid
    % The rock (perm/poro/ntg)
    rock

    % Maximum relative pressure change
    dpMaxRel
    % Maximum absolute pressure change
    dpMaxAbs
    % Maximum Relative saturation change
    dsMaxRel
    % Maximum absolute saturation change
    dsMaxAbs

    % Maximum pressure allowed in reservoir
    maximumPressure
    % Minimum pressure allowed in reservoir
    minimumPressure

    % Indicator showing if the aqua/water phase is present
    water
    % Indicator showing if the gas phase is present
    gas
    % Indicator showing if the oil phase is present
    oil
    % Names of primary variables interpreted as saturations, i.e. so
    % that they will sum to one when updated.
    saturationVarNames
    % Names of components
    componentVarNames
    
    % Use alternate tolerance scheme
    useCNVConvergence

    % CNV tolerance (inf-norm-like)
    toleranceCNV;

    % MB tolerance values (2-norm-like)
    toleranceMB;
    % Well tolerance if CNV is being used
    toleranceWellBHP;
    % Well tolerance if CNV is being used
    toleranceWellRate;

    % Input data used to instantiate the model
    inputdata
    % Add extra output to wellsol/states for derived quantities
    extraStateOutput
    % Output extra information for wells
    extraWellSolOutput
    % Output fluxes
    outputFluxes
    % Upstream weighting of injection cells
    upstreamWeightInjectors
    
    % Vector for the gravitational force
    gravity
    
    % Well model used to compute fluxes etc from well controls
    wellmodel
end

methods
    % --------------------------------------------------------------------%
    function model = ReservoirModel(G, varargin)
        model = model@PhysicalModel(G);

        if nargin == 1 || ischar(varargin{1})
            % We were given only grid + any keyword arguments
            doSetup = false;
        else
            assert(nargin >= 3)
            % We are being called in format
            % ReservoirModel(G, rock, fluid, ...)
            model.rock  = varargin{1};
            model.fluid = varargin{2};
            
            % Rest of arguments should be keyword/value pairs.
            varargin = varargin(3:end);
            % We have been provided the means, so we will execute setup
            % phase after parsing other inputs and defaults.
            doSetup = true;
        end
        
        model.dpMaxRel = inf;
        model.dpMaxAbs = inf;

        model.minimumPressure = -inf;
        model.maximumPressure =  inf;

        model.dsMaxAbs = .2;
        model.dsMaxRel = inf;

        model.nonlinearTolerance = 1e-6;
        model.inputdata = [];

        model.useCNVConvergence = false;
        model.toleranceCNV = 1e-3;
        model.toleranceMB = 1e-7;
        model.toleranceWellBHP = 1*barsa;
        model.toleranceWellRate = 1/day;

        model.saturationVarNames = {'sw', 'so', 'sg'};
        model.componentVarNames  = {};

        model.extraStateOutput = false;
        model.extraWellSolOutput = true;
        model.outputFluxes = true;
        model.upstreamWeightInjectors = false;
        % Gravity defaults to the global variable
        model.gravity = gravity(); %#ok
        [model, unparsed] = merge_options(model, varargin{:}); %#ok

        % Base class does not support any phases
        model.water = false;
        model.gas = false;
        model.oil = false;
        
        if doSetup
            if isempty(G) || isempty(model.rock)
                dispif(model.verbose, 'mrst:ReservoirModel', ...
                    'Invalid grid/rock pair supplied. Operators have not been set up.')
            else
                model.operators = setupOperatorsTPFA(G, model.rock, 'deck', model.inputdata);
            end
        end
    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@PhysicalModel(model, state);
        active = model.getActivePhases();
        nPh = nnz(active);
        nc = model.G.cells.num;
        model.checkProperty(state, 'Pressure', [nc, 1], [1, 2]);
        if nPh > 1
            model.checkProperty(state, 'Saturation', [nc, nPh], [1, 2]);
        end
        state = model.wellmodel.validateState(state);
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        % Generic update function for reservoir models containing wells

        % Split variables into three categories: Regular/rest variables, saturation
        % variables (which sum to 1 after updates) and well variables (which live
        % in wellSol and are in general more messy to work with).
        [restVars, satVars] = model.splitPrimaryVariables(problem.primaryVariables);

        % Update the wells
        if isfield(state, 'wellSol')
            [state.wellSol, restVars] = model.wellmodel.updateWellSol(state.wellSol, problem, dx, drivingForces, restVars);
        end
        
        % Update saturations in one go
        state  = model.updateSaturations(state, dx, problem, satVars);

        if ~isempty(restVars)
            % Handle pressure seperately
            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
            restVars = model.stripVars(restVars, 'pressure');

            % Update remaining variables (tracers, temperature etc)
            for i = 1:numel(restVars);
                 state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
            end
        end

        report = [];
    end

    % --------------------------------------------------------------------%
    function model = setupOperators(model, G, rock, varargin)
        % Set up divergence/gradient/transmissibility operators
        model.operators = setupOperatorsTPFA(G, rock, varargin{:});
    end

    % --------------------------------------------------------------------%
    function [convergence, values, names] = checkConvergence(model, problem, varargin)
        if model.useCNVConvergence
            % Use convergence model similar to commercial simulator
            [conv_cells, v_cells, isWOG, namesWOG] = CNV_MBConvergence(model, problem);
            [conv_wells, v_wells, isWell, namesWell] = checkWellConvergence(model, problem);

            % Get the values for all equations, just in case there are some
            % values that are not either wells or standard 3ph conservation
            % equations
            values_all = norm(problem, inf);
            rest = ~(isWOG | isWell);
            
            tol = model.nonlinearTolerance;
            convergence = all(conv_cells) && ...
                          all(conv_wells) && ...
                          all(values_all(rest) < tol);
                      
            values = [v_cells, v_wells];
            names = horzcat(namesWOG, namesWell);
        else
            % Use strict tolerances on the residual without any 
            % fingerspitzengefuhlen by calling the parent class
            [convergence, values, names] = checkConvergence@PhysicalModel(model, problem, varargin{:});
        end            
    end

    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
        switch(lower(name))
            case {'t', 'temperature'}
                fn = 'T';
                index = 1;
            case {'sw', 'water'}
                index = model.satVarIndex('sw');
                fn = 's';
            case {'so', 'oil'}
                index = model.satVarIndex('so');
                fn = 's';
            case {'sg', 'gas'}
                index = model.satVarIndex('sg');
                fn = 's';
            case {'s', 'sat', 'saturation'}
                index = ':';
                fn = 's';
            case {'pressure', 'p'}
                index = 1;
                fn = 'pressure';
            case 'wellsol'
                % Use colon to get all variables, since the wellsol may
                % be empty
                index = ':';
                fn = 'wellSol';
            otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@PhysicalModel(model, name);
        end
    end
    
    % --------------------------------------------------------------------%
    function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@PhysicalModel(model);
        % Support for wells
        forces.W   = [];
        % Support for boundary conditions
        forces.bc  = [];
        % Support for direct source terms
        forces.src = [];
    end

    % --------------------------------------------------------------------%
    function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
        % Split a set of primary variables into three groups:
        % Well variables, saturation variables and the rest. This is
        % useful because the saturation variables usually are updated
        % together, and the well variables are a special case.
        wellvars = model.wellmodel.getPrimaryVariableNames();
        isSat   = cellfun(@(x) any(strcmpi(model.saturationVarNames, x)), vars);
        isWells = cellfun(@(x) any(strcmpi(wellvars, x)), vars);

        wellVars = vars(isWells);
        satVars  = vars(isSat);

        restVars = vars(~isSat & ~isWells);
    end

    % --------------------------------------------------------------------%
    function [isActive, phInd] = getActivePhases(model)
        % Get active flag for the canonical phase ordering (water, oil
        % gas as on/off flags).
        isActive = [model.water, model.oil, model.gas];
        if nargout > 1
            phInd = find(isActive);
        end
    end

    % --------------------------------------------------------------------%
    function phNames = getPhaseNames(model)
        % Get the active phases in canonical ordering
        tmp = 'WOG';
        active = model.getActivePhases();
        phNames = tmp(active);
    end

    % --------------------------------------------------------------------%
    function i = getPhaseIndex(model, phasename)
        % Query the index of a phase ('W', 'O', 'G')
        active = model.getPhaseNames();
        i = find(active == phasename);
    end 

    % --------------------------------------------------------------------%
    function state = updateSaturations(model, state, dx, problem, satVars)
        % Update saturations (likely state.s) under the constraint that
        % the sum of volume fractions is always equal to 1. This
        % assumes that we have solved for n - 1 phases when n phases
        % are present.
        if nargin < 5
            % Get the saturation names directly from the problem
            [~, satVars] = ...
                splitPrimaryVariables(model, problem.primaryVariables);
        end
        if isempty(satVars)
            % No saturations passed, nothing to do here.
            return
        end
        % Solution variables should be saturations directly, find the missing
        % link
        saturations = lower(model.saturationVarNames);
        fillsat = setdiff(saturations, lower(satVars));
        assert(numel(fillsat) == 1)
        fillsat = fillsat{1};

        % Fill component is whichever saturation is assumed to fill up the rest of
        % the pores. This is done by setting that increment equal to the
        % negation of all others so that sum(s) == 0 at end of update
        solvedFor = ~strcmpi(saturations, fillsat);
        ds = zeros(model.G.cells.num, numel(saturations));

        tmp = 0;
        for i = 1:numel(saturations)
            if solvedFor(i)
                v = model.getIncrement(dx, problem, saturations{i});
                ds(:, i) = v;
                % Saturations added for active variables must be subtracted
                % from the last phase
                tmp = tmp - v;
            end
        end
        ds(:, ~solvedFor) = tmp;
        % We update all saturations simultanously, since this does not bias the
        % increment towards one phase in particular.
        state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);

        % Ensure that values are within zero->one interval, and
        % re-normalize if any values were capped
        bad = any((state.s > 1) | (state.s < 0), 2);
        if any(bad)
            state.s(bad, :) = min(state.s(bad, :), 1);
            state.s(bad, :) = max(state.s(bad, :), 0);
            state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), sum(state.s(bad, :), 2));
        end
    end

    % --------------------------------------------------------------------%
    function state = setPhaseData(model, state, data, fld, subs)
        % Given a structure field name and a cell array of data for
        % each phase, store them as columns with the given field name
        if nargin == 4
            subs = ':';
        end
        isActive = model.getActivePhases();

        ind = 1;
        for i = 1:numel(data)
            if isActive(i)
                state.(fld)(subs, ind) = data{i};
                ind = ind + 1;
            end
        end
    end
    
    % --------------------------------------------------------------------%
    function state = storeFluxes(model, state, vW, vO, vG)
        % Utility function for storing the interface fluxes in the state
        isActive = model.getActivePhases();

        internal = model.operators.internalConn;
        state.flux = zeros(numel(internal), sum(isActive));
        phasefluxes = {double(vW), double(vO), double(vG)};
        state = model.setPhaseData(state, phasefluxes, 'flux', internal);
    end
    
    function state = storeBoundaryFluxes(model, state, qW, qO, qG, forces)
        if isempty(forces.bc) || ~isfield(forces, 'bc')
            return
        end
        phasefluxes = {double(qW), double(qO), double(qG)};
        state = model.setPhaseData(state, phasefluxes, 'flux', forces.bc.face);
    end
    
    % --------------------------------------------------------------------%
    function state = storeMobilities(model, state, mobW, mobO, mobG)
        % Utility function for storing the mobilities in the state
        isActive = model.getActivePhases();

        state.mob = zeros(model.G.cells.num, sum(isActive));
        mob = {double(mobW), double(mobO), double(mobG)};
        state = model.setPhaseData(state, mob, 'mob');
    end
    
    % --------------------------------------------------------------------%
    function state = storeUpstreamIndices(model, state, upcw, upco, upcg)
        % Store upstream indices, so that they can be reused for other
        % purposes.
        isActive = model.getActivePhases();

        nInterfaces = size(model.operators.N, 1);
        state.upstreamFlag = false(nInterfaces, sum(isActive));
        mob = {upcw, upco, upcg};
        state = model.setPhaseData(state, mob, 'upstreamFlag');
    end
    
    % --------------------------------------------------------------------%
    function state = storeDensity(model, state, rhoW, rhoO, rhoG)
        % Store compressibility / surface factors for plotting and
        % output.
        isActive = model.getActivePhases();

        state.rho = zeros(model.G.cells.num, sum(isActive));
        rho = {double(rhoW), double(rhoO), double(rhoG)};
        state = model.setPhaseData(state, rho, 'rho');
    end
    % --------------------------------------------------------------------%
    function state = storebfactors(model, state, bW, bO, bG)
        % Store compressibility / surface factors for plotting and
        % output.
        isActive = model.getActivePhases();

        state.bfactor = zeros(model.G.cells.num, sum(isActive));
        b = {double(bW), double(bO), double(bG)};
        state = model.setPhaseData(state, b, 'bfactor');
    end
    
    % --------------------------------------------------------------------%
    function i = satVarIndex(model, name)
        % Find the index of a saturation variable by name
        i = find(strcmpi(model.saturationVarNames, name));
    end
    
    % --------------------------------------------------------------------%
    function i = compVarIndex(model, name)
        % Find the index of a component variable by name
        i = find(strcmpi(model.componentVarNames, name));
    end
    
    % --------------------------------------------------------------------%
    function varargout = evaluateRelPerm(model, sat, varargin)
        % Evaluate the fluid relperm. Depending on the active phases,
        % we must evaluate the right fluid relperm functions and
        % combine the results. This function calls the appropriate
        % static functions.
        active = model.getActivePhases();
        nph = sum(active);
        assert(nph == numel(sat), ...
        'The number of saturations must equal the number of active phases.')
        varargout = cell(1, nph);
        names = model.getPhaseNames();

        if nph > 1
            fn = ['relPerm', names];
            [varargout{:}] = model.(fn)(sat{:}, model.fluid, varargin{:});
        elseif nph == 1
            % Call fluid interface directly if single phase
            varargout{1} = model.fluid.(['kr', names])(sat{:}, varargin{:});
        end
    end
    
    % --------------------------------------------------------------------%
    function g = getGravityVector(model)
        % Get the gravity vector used to instantiate the model
        if isfield(model.G, 'griddim')
            dims = 1:model.G.griddim;
        else
            dims = ':';
        end
        g = model.gravity(dims);
    end
    
    % --------------------------------------------------------------------%
    function gdxyz = getGravityGradient(model)
        % Get gradient of gravity on the faces
        assert(isfield(model.G, 'cells'), 'Missing cell field on grid');
        assert(isfield(model.G.cells, 'centroids'),...
            'Missing centroids field on grid. Consider using computeGeometry first.');
        
        g = model.getGravityVector();
        gdxyz = model.operators.Grad(model.G.cells.centroids) * g';
    end
    
% --------------------------------------------------------------------%
    function scaling = getScalingFactorsCPR(model, problem, names)%#ok
        % Return cell array of scaling factors for approximate pressure
        % equation in CPR preconditioner.
        %
        % Either one value per equation, or one cell wise value per
        % equation.
        %
        % Scaling should be the size of the names sent in.
        scaling = cell(numel(names), 1);
        [scaling{:}] = deal(1);
    end

end

methods (Static)
    % --------------------------------------------------------------------%
    function [krW, krO, krG] = relPermWOG(sw, so, sg, f, varargin)
        % Three phase, water / oil / gas relperm.
        swcon = 0;
        if isfield(f, 'sWcon')
            swcon = f.sWcon;
        end
        swcon = min(swcon, double(sw)-1e-5);

        d  = (sg+sw-swcon);
        ww = (sw-swcon)./d;
        krW = f.krW(sw, varargin{:});

        wg = 1-ww;
        krG = f.krG(sg, varargin{:});

        krow = f.krOW(so, varargin{:});
        krog = f.krOG(so,  varargin{:});
        krO  = wg.*krog + ww.*krow;
    end
    
    % --------------------------------------------------------------------%
    function [krW, krO] = relPermWO(sw, so, f, varargin)
        % Two phase oil-water relperm
        krW = f.krW(sw, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOW(so, varargin{:});
        end
    end
    
    % --------------------------------------------------------------------%
    function [krO, krG] = relPermOG(so, sg, f, varargin)
        % Two phase oil-gas relperm.
        krG = f.krG(sg, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOG(so, varargin{:});
        end
    end
    
    % --------------------------------------------------------------------%
    function [krW, krG] = relPermWG(sw, sg, f, varargin)
        % Two phase water-gas relperm
        krG = f.krG(sg, varargin{:});
        krW = f.krW(sw, varargin{:});
    end
    
    % --------------------------------------------------------------------%
    function ds = adjustStepFromSatBounds(s, ds)
        % Ensure that cellwise increment for each phase is done with
        % the same length, in a manner that avoids saturation
        % violations.
        tmp = s + ds;

        violateUpper =     max(tmp - 1, 0);
        violateLower = abs(min(tmp    , 0));

        violate = max(violateUpper, violateLower);

        [worst, jj]= max(violate, [], 2);

        bad = worst > 0;
        if any(bad)
            w = ones(size(s, 1), 1);
            for i = 1:size(s, 2)
                ind = bad & jj == i;
                dworst = abs(ds(ind, i));

                w(ind) = (dworst - worst(ind))./dworst;
            end
            ds(bad, :) = bsxfun(@times, ds(bad, :), w(bad, :));
        end
    end
end
end

