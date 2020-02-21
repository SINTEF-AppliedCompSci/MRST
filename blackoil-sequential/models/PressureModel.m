classdef PressureModel < WrapperModel
    properties
        incTolPressure = 1e-3;
        useIncTol = true;
        reductionStrategy;
        pressureIncTolType = 'relative';
        pressureTol = inf;
    end
    
    methods
        function model = PressureModel(parent, varargin)
            if isprop(parent, 'useCNVConvergence')
                parent.useCNVConvergence = false;
            end
            model = model@WrapperModel(parent);
            model = merge_options(model, varargin{:});
            if isempty(model.reductionStrategy)
                if isa(model, 'ThreePhaseBlackOilModel')
                    s = 'analytic';
                else
                    s = 'numerical';
                end
                model.reductionStrategy = s;
            end
            model.AutoDiffBackend = parent.AutoDiffBackend;
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            % Get the AD state for this model
            [vars, names, origin] = model.getPrimaryVariables(state);
            switch model.reductionStrategy
                case {'analytic', 'numerical'}
                    isP = strcmp(names, 'pressure');
                    origP = origin{isP};
                    keep = isP | cellfun(@(x) ~strcmp(x, origP), origin);
                otherwise
                    error('Unknown reduction strategy %s', model.reductionStrategy);
            end
            if init
                [vars{keep}] = model.AutoDiffBackend.initVariablesAD(vars{keep});
            end
            state = model.initStateAD(state, vars, names, origin);
            % Not all were AD-initialized
            names = names(keep);
            origin = origin(keep);
            assert(strcmpi(names{1}, 'pressure'));
        end

        function [eqs, names, types, state] = getModelEquations(pmodel, state0, state, dt, drivingForces)
            model = pmodel.parentModel;
            ncomp = model.getNumberOfComponents();
            [eqs, names, types, state] = model.getModelEquations(state0, state, dt, drivingForces);
            ceqs = eqs(1:ncomp);
            
            w = model.getProp(state, 'PressureReductionFactors');
            % Assemble equations and add in sources
            pressure_equation = 0;
            for i = 1:numel(ceqs)
                pressure_equation = pressure_equation + w{i}.*ceqs{i};
            end
            subs = (ncomp+1):numel(eqs);
            eqs = [{pressure_equation}, eqs(subs)];
            names = ['pressure', names(subs)];
            types = ['cell', types(subs)];
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@WrapperModel(model, state, state0, dt, drivingForces);
            state = assignReductionFactorProps(model, state, state0, dt);
            % Ensure that saturations are normalized
            state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            if strcmpi(model.reductionStrategy, 'numerical')
                state = updateReductionFactorProps(model, state);
            end
            p0 = state.pressure;
            [state, report] = model.parentModel.updateState(state, problem, dx, drivingForces);
            switch lower(model.pressureIncTolType)
                case 'relative'
                    range = max(p0) - min(p0);
                    if range == 0
                        range = 1;
                    end
                    dp = (state.pressure - p0)./range;
                case 'absolute'
                    dp = state.pressure - p0;
                case 'norm'
                    dp = (state.pressure - p0)/norm(state.pressure, inf);
                otherwise
                    error('Unknown pressure increment %s', model.pressureIncTolType);
            end
            state.pressureChange = dp;
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
            ptol = model.pressureTol;
            itol = model.incTolPressure;
            usePTol = isfinite(ptol);
            useITol = isfinite(itol) && model.useIncTol;
            assert(usePTol || useITol, 'Pressure model has no valid way of checking convergence!');
            if usePTol
                % Check the actual value of the pressure equation. Note:
                % Requires that the equation is actually properly scaled!
                % The scaling is implicit from the pressure reduction
                % factors.
                convergence(1) = values(1) < ptol;
            else
                % Skip first equation
                names = names(2:end);
                convergence = convergence(2:end);
                values = values(2:end);
            end
            
            if useITol
                % Check the increment tolerance  for pressure
                if problem.iterationNo > 1 && ~isnan(problem.iterationNo)
                    dp = norm(problem.state.pressureChange, inf);
                else
                    dp = inf;
                end
                convergence = [dp < itol, convergence];
                names = ['Delta P', names];
                values = [dp, values];
            end
        end

        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@WrapperModel(model, varargin{:});
            % Clean up internal fields
            if isfield(state, 'statePressure')
                state = rmfield(state, 'statePressure');
            end
            if isfield(state, 'sT')
                state.sT = sum(value(state.s), 2);
            end
            if isfield(state, 'pressureChange')
                state = rmfield(state, 'pressureChange');
            end
            if isfield(state, 'reductionFactorProps')
                state = rmfield(state, 'reductionFactorProps');
            end 
            state.statePressure = state;
        end
        
        function rhoS = getSurfaceDensities(model)
            rhoS = model.parentModel.getSurfaceDensities();
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            pmodel = model.parentModel;
            pvt = pmodel.PVTPropertyFunctions;
            switch model.reductionStrategy
                case 'analytic'
                    assert(isa(pmodel, 'ThreePhaseBlackOilModel'), ...
                        'Analytical pressure reduction factors currently only implemented for black-oil.');
                case 'numerical'
                    pvt = pvt.setStateFunction('PressureReductionFactors', NumericalPressureReductionFactors(pmodel));
                otherwise
                    error('Unknown reduction strategy');
            end
            model.parentModel.PVTPropertyFunctions = pvt;
        end
        
    end
end

function state = assignReductionFactorProps(model, state, state0, dt)
    mass0 = model.parentModel.getProp(state0, 'ComponentTotalMass');
    props = struct('mass0'   , {mass0}, ...
                   'dt'      , dt     , ...
                   'pressure', []     , ...
                   'weights' , []     );
    state.reductionFactorProps = props;
end

function state = updateReductionFactorProps(model, state)
    pressure = model.getProp(state, 'pressure');
    weights  = model.parentModel.getProp(state, 'PressureReductionFactors');
    state.reductionFactorProps.pressure = pressure;
    state.reductionFactorProps.weights = horzcat(weights{:});
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
