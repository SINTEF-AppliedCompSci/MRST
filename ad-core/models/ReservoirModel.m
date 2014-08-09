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
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
        
        % Maximum relative pressure change
        dpMax
        % Maximum absolute saturation change
        dsMax
        
        % Water phase present
        water
        % Gas phase present
        gas
        % Oil phase present
        oil

        
        % Names of primary variables interpreted as saturations, i.e. so
        % that they will sum to one when updated.
        saturationVarNames
        % Names of well fields that may be updated by the model.
        wellVarNames
        
        % Input data used to instantiate the model
        inputdata
    end
    
    methods
        function model = ReservoirModel(G, rock, fluid, varargin) %#ok
            model = model@PhysicalModel(G);
            
            model.dpMax = inf;
            model.dsMax = .2;
            model.nonlinearTolerance = 1e-6;
            model.inputdata = [];
            
            model.saturationVarNames = {'sw', 'so', 'sg'};
            model.wellVarNames = {'qWs', 'qOs', 'qGs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;
            
            % Physical model
            model.fluid = fluid;
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Generic update function for reservoir models containing wells

            % Split variables into three categories: Regular/rest variables, saturation
            % variables (which sum to 1 after updates) and well variables (which live
            % in wellSol and are in general more messy to work with).
            [restVars, satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);
            
            % Update saturations in one go
            state  = model.updateSaturations(state, dx, problem, satVars);
            
            if ~isempty(restVars)
                % Handle pressure seperately
                state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMax);
                restVars = model.stripVars(restVars, 'pressure');

                % Update remaining variables (tracers, temperature etc)
                for i = 1:numel(restVars);
                     state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
                end
            end

            % Update the wells
            state.wellSol = model.updateWellSol(state.wellSol, dx, problem, wellVars);
            report = [];
        end

        function model = setupOperators(model, G, rock, varargin)
            % Set up divergence/gradient/transmissibility operators
            if isempty(G) || isempty(rock)
                warning('mrst:ReservoirModel', ...
                'Invalid grid/rock pair supplied. Operators have not been set up.')
                return;
            end
            model.operators = setupOperatorsTPFA(G, rock, varargin{:});
        end
                   
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            % Setup and pass on driving forces
            vararg = {};
            driving = struct('Wells', [], 'bc', [], 'src', []);
            
            if isfield(control, 'W') && ~isempty(control.W)
                vararg = [vararg, 'Wells', control.W];
                driving.Wells = control.W;
            end

            if isfield(control, 'bc') && ~isempty(control.bc)
                vararg = [vararg, 'bc', control.bc];
                driving.bc = control.bc;
            end
            
            if isfield(control, 'src') && ~isempty(control.src)
                vararg = [vararg, 'src', control.src];
                driving.src = control.src;
            end
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'t', 'temperature'}
                    fn = 'T';
                    index = 1;
                case {'sw', 'water'}
                    index = find(strcmpi(model.saturationVarNames, 'sw'));
                    fn = 's';
                case {'so', 'oil'}
                    index = find(strcmpi(model.saturationVarNames, 'so'));
                    fn = 's';
                case {'sg', 'gas'}
                    index = find(strcmpi(model.saturationVarNames, 'sg'));
                    fn = 's';
                case {'s', 'sat', 'saturation'}
                    index = 1:numel(model.saturationVarNames);
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
        
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            isSat   = cellfun(@(x) any(strcmpi(model.saturationVarNames, x)), vars);
            isWells = cellfun(@(x) any(strcmpi(model.wellVarNames, x)), vars);
            
            wellVars = vars(isWells);
            satVars  = vars(isSat);
                
            restVars = vars(~isSat & ~isWells);
        end
        
        function [isActive, phInd] = getActivePhases(model)
            isActive = [model.water, model.oil, model.gas];
            if nargout > 1
                phInd = find(isActive);
            end
        end
        
        function phNames = getPhaseNames(model)
            tmp = 'WOG';
            active = model.getActivePhases();
            phNames = tmp(active);
        end
        
        function i = getPhaseIndex(model, phasename)
            active = model.getPhaseNames();
            i = find(active == phasename);
        end 
        
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
            state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMax);
        end
        
        function wellSol = updateWellSol(model, wellSol, dx, problem, wellVars)
            % Update the wellSol struct
            if nargin < 4
                % Get the well variables directly from the problem,
                % otherwise assume that they are known by the user
                [~, ~, wellVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            
            for i = 1:numel(wellVars)
                wf = wellVars{i};
                dv = model.getIncrement(dx, problem, wf);

                if strcmpi(wf, 'bhp')
                    % Bottom hole is a bit special - we apply the pressure update
                    % limits here as well.
                    bhp = vertcat(wellSol.bhp);
                    dv = model.limitUpdateRelative(dv, bhp, model.dpMax);
                end

                for j = 1:numel(wellSol)
                    wellSol(j).(wf) = wellSol(j).(wf) + dv(j);
                end
            end
        end
    end

    methods (Static)

    end

end

