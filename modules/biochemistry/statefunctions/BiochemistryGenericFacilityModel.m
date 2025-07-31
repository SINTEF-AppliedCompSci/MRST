classdef BiochemistryGenericFacilityModel < GenericFacilityModel
    % Generic facility model for biochemistry simulations
    
    properties
        bacterialFormulation = 'bacterialmodel';  % Formulation for bacterial transport
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, useDefaults)
            % Set up state function groupings using parent, and add
            % model-specific functionality
            
            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end
            
            % Set up state function groupings using parent
            model = setupStateFunctionGroupings@GenericFacilityModel(model, useDefaults);
            
            % Add biochemistry-specific state functions
            ffd = model.FacilityFlowDiscretization;
            ffd = ffd.setStateFunction('BacterialMass', BacterialMass(model));
            ffd = ffd.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
            ffd = ffd.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
            ffd = ffd.setStateFunction('BactConvRate', BactConvertionRate(model));
            model.FacilityFlowDiscretization = ffd;
        end
        
        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            % Initialize AD state from double state
            state = initStateAD@GenericFacilityModel(model, state, vars, names, origin);
        end
        
        %-----------------------------------------------------------------%
        function names = getBasicPrimaryVariableNames(model)
            % Get names of primary variables
            
            % Get parent class primary variable names
            names = getBasicPrimaryVariableNames@GenericFacilityModel(model);
            
            % Check if we have bacterial primary variables
            if strcmpi(model.bacterialFormulation, 'none') || ...
               strcmpi(model.primaryVariableSet, 'none')
                return
            end
        end
        
        %-----------------------------------------------------------------%
        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            % Get primary variables
            [variables, names, map] = ...
                getBasicPrimaryVariables@GenericFacilityModel(model, wellSol);
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get field name for variable
            [fn, index] = getVariableField@GenericFacilityModel(model, name, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function src_growthdecay = getBacteriaSources(model, fd, state, state0, dt)
            % Compute bacterial growth and decay sources
            
            flowState = fd.buildFlowState(model, state, state0, dt);  
            psigrowth = model.getProps(flowState, 'PsiGrowthRate');
            psidecay = model.getProps(flowState, 'PsiDecayRate');
           % bmass = model.getProps(flowState, 'BacterialMass');

            src_growthdecay = (psigrowth - psidecay);
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Get model equations including parent equations
            [eqs, names, types, state] = ...
                getModelEquations@GenericFacilityModel(model, state0, state, dt, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function [values, tolerances, names, evaluated] = getFacilityConvergenceValues(model, problem, varargin)
            % Get convergence values for facility
            [values, tolerances, names, evaluated] = ...
                getFacilityConvergenceValues@GenericFacilityModel(model, problem, varargin{:});
        end
    end
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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