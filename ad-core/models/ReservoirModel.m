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
        
        % Names of each component
        componentNames
        
        % Names of each saturation variables, corresponding to their order in state.s
        saturationNames
        
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
            
            model.saturationNames = {'sw', 'so', 'sg'};
            model.componentNames = {};
            
            model = merge_options(model, varargin{:});
            
            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;
            
            % Physical model
            model.fluid = fluid;
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
            index = [];
            switch(lower(name))
                case {'t', 'temperature'}
                    fn = 'T';
                case {'sw', 'water'}
                    index = find(strcmpi(model.saturationNames, 'sw'));
                    fn = 's';
                case {'so', 'oil'}
                    index = find(strcmpi(model.saturationNames, 'so'));
                    fn = 's';
                case {'sg', 'gas'}
                    index = find(strcmpi(model.saturationNames, 'sg'));
                    fn = 's';
                case {'s', 'sat', 'saturation'}
                    index = 1:numel(model.saturationNames);
                    fn = 's';
                case {'pressure', 'p'}
                    index = 1;
                    fn = 'pressure';
                case 'wellsol'
                    index = 1;
                    fn = 'wellSol';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name);
            end
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
    end

    methods (Static)

    end

end

