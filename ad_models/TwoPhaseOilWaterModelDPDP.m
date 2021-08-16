classdef TwoPhaseOilWaterModelDPDP < TwoPhaseOilWaterModelDP
    % Dual porosity - dual permeability two phase oil/water model
    properties
        % the matrix operators
        operators_matrix
    end

    methods
        function model = TwoPhaseOilWaterModelDPDP(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModelDP(G, rock, fluid, varargin{:});

            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;

            % Blackoil -> use CNV style convergence
            model.useCNVConvergence = true;

            % Comment for 2020a
            %model.matrixVarNames = {'pom', 'swm'};

            model.extraStateOutput = 1;

            % model = merge_options(model, varargin{:});
            
            % Set up the operators for the matrix
            model.operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, 'deck', model.inputdata);
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterDPDP(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});

        end
        
%         % --------------------------------------------------------------------%
%         function rhoS = getSurfaceDensities(model)
%             % Get the surface densities of the active phases in canonical
%             % ordering (WOG, with any inactive phases removed).
%             %
%             % RETURNS:
%             %   rhoS - 1 x n double array of surface densities.
%             names = model.getPhaseNames();
%             rhoS = arrayfun(@(x) model.fluid.(['rho', x, 'S']), names);
%             rhoS = [rhoS, rhoS];    % Quick fix
%         end
               
        
%         %--------------------------------------------------------------------%
%         function [phNames, longNames] = getPhaseNames(model)
%             Get short and long names of the present phases.
%             
%             SYNOPSIS:
%               [phNames, longNames] = model.getPhaseNames();
%             
%             PARAMETERS:
%               model    - Class instance
%             
%             RETURNS:
%               phNames   - Cell array containing the short hanes ('W', 'O',
%                           G') of the phases present
%               
%               longNames - Longer names ('water', 'oil', 'gas') of the phases
%                           present.            
%                        
%             phNames = 'WO';
%             longNames = {'water', 'oil', 'water_matrix', 'oil_matrix'};
%         end


%         %--------------------------------------------------------------------%      
%         function names = getComponentNames(model)
%             names = getComponentNames@TwoPhaseOilWaterDPModel(model);
%             names{end+1} = 'water_matrix';
%             names{end+1} = 'oil_matrix';
%         end
        
        
        %--------------------------------------------------------------------%  
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            % For a given component conservation equation, compute and add in
            % source terms for a specific source/bc where the fluxes have
            % already been computed.
            %
            % PARAMETERS:
            %
            %   model  - (Base class, automatic)
            %
            %   cname  - Name of the component. Must be a property known to the
            %            model itself through `getProp` and `getVariableField`.
            %
            %   eq     - Equation where the source terms are to be added. Should
            %            be one value per cell in the simulation grid (model.G)
            %            so that the src.sourceCells is meaningful.
            %
            %   component - Cell-wise values of the component in question. Used
            %               for outflow source terms only.
            %
            %   src    - Source struct containing fields for fluxes etc. Should
            %            be constructed from force and the current reservoir
            %            state by `computeSourcesAndBoundaryConditionsAD`.
            %
            %   force  - Force struct used to produce src. Should contain the
            %            field defining the component in question, so that the
            %            inflow of the component through the boundary condition
            %            or source terms can accurately by estimated.
            
            
%             % test
%             cells = src.sourceCells;
%             qW = src.phaseMass{1}./model.fluid.rhoWS;
%             qW.val = ones( length(qW.val), 1);
%             qC = qW;
%             %isInj = qW > 0;
%             %qC = (isInj.*1 + ~isInj.*1).*qW;
%             eq(cells) = eq(cells) - qC;            
%             src.components{end+1} = qC; 
            
            %{
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'polymer'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;                    
            %}
        end


    end
    
end

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
 