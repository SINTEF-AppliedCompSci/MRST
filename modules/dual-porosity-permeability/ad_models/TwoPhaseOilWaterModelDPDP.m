classdef TwoPhaseOilWaterModelDPDP < ThreePhaseBlackOilModelDP
    % The dual porosity - dual permeability two phase oil/water model
    properties
        % the matrix operators
        operators_matrix
    end

    methods
        function model = TwoPhaseOilWaterModelDPDP(G, rock, fluid, varargin)
            % NA: pass varargin to the parent constructor e.g. for eventual
            % transmissibility multipliers            
            model = model@ThreePhaseBlackOilModelDP(G, rock, fluid, varargin{:});

            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;

            % Blackoil -> use CNV style convergence
            model.useCNVConvergence = true;

            model.extraStateOutput = 1;
            
            % Set up the matrix and fracture operators with potentially 
            % provided options
            if isempty(model.inputdata)
                model.operators = setupOperatorsTPFA(G, model.rock);
                model.operators_matrix = setupOperatorsTPFA(G, model.rock_matrix); 
                operators_matrix = setupOperatorsTPFA(G, model.rock_matrix); 
            else
                % Interpret inputdata as an Eclipse deck structure with
                % the additional fields 'neighbors' 'trans', and 'porv'.
                fields = {'neighbors', 'trans', 'porv'};
                opts = {};
                opts_matrix = {};
                for i = 1:numel(fields)
                    if isfield(model.inputdata, fields{i})
                        opts{end + 1} = fields{i};
                        opts{end + 1} = model.inputdata.(fields{i}){1};
                        opts_matrix{end + 1} = fields{i};
                        opts_matrix{end + 1} = model.inputdata.(fields{i}){2};
                    end
                end
                model.operators = setupOperatorsTPFA(G, model.rock, opts{:});
                operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, opts_matrix{:});    
                model.operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, opts_matrix{:}); 
            end
            model.operators.pv_matrix = operators_matrix.pv;            
   
            model = merge_options(model, varargin{:});

        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterDPDP(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});

        end
             
        
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
           
        end
        
       
        % --------------------------------------------------------------------%
        function varargout = evaluateRelPerm(model, sat, varargin)
            % Overwrite the parent class' method
            active = model.getActivePhases();
            nph = sum(active);
            assert(nph == numel(sat), ...
            'The number of saturations must equal the number of active phases.')
            varargout = cell(1, nph);
            names = model.getPhaseNames();

            % Ensure that both fracture and matrix relative permeabilities are evaluated
            if nargin == 4
                if strcmpi(varargin{1}, 'medium') && strcmpi(varargin{2}, 'fracture')
                    f = model.fluid;
                elseif strcmpi(varargin{1}, 'medium') && strcmpi(varargin{2}, 'matrix')
                    f = model.fluid_matrix;
                else
                    error(['Unsupported options ' varargin{1} ', ' varargin{2} ...
                           ' in evaluateRelPerm']);
                end
            else
                error('Medium not specified in evaluateRelPerm');
            end

            if nph > 1
                fn = ['relPerm', names];
                [varargout{:}] = model.(fn)(sat{:}, f, varargin{:});
            elseif nph == 1
                % Call fluid interface directly if single phase
                varargout{1} = model.fluid.(['kr', names])(sat{:}, varargin{:});
            end
        end        
        
        
        %--------------------------------------------------------------------%  
        function state = validateState(model, state)
            % Dummy validation to be called in simulateScheduleAD.
            return
        end        


    end
    
end

%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

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
 