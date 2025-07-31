classdef DiffusiveBactFlux < StateFunction & UpwindProperty
    % Bacterial diffusive flux computation
    %
    % SYNOPSIS:
    %   df = DiffusiveBactFlux(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the diffusive flux of bacterial cells in the liquid phase,
    %   accounting for:
    %   - Bacterial concentration gradients
    %   - Liquid phase saturation
    %   - Phase densities
    %   - Effective diffusion coefficient
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with bacterial modeling enabled
    %
    % OPTIONAL PARAMETERS:
    %   'Db' - Diffusion coefficient (default: 1e-8 m^2/s)
    %   'upwind_name' - Name of upwind flag function (default: 'PhaseUpwindFlag')
    %
    % RETURNS:
    %   Class instance ready for use in simulation
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional

    properties (Access = protected)
        Db = 1e-8*meter^2/second; % Effective diffusion coefficient [m^2/s]
        upwind_name; % Name of state function providing upwind flag
    end
    
    methods
        function gp = DiffusiveBactFlux(model, varargin)                                  
            % Constructor for bacterial diffusive flux
            
            % Process optional inputs
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            
            % Initialize parent classes
            gp@StateFunction(model);
            gp@UpwindProperty(upstr);
            
            % Set dependencies
            gp = gp.dependsOn(gp.upwind_name);
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'Density'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'s'}, 'state');
            
            gp.label = '\mathbf{F}_{bio}^{diff}'; % LaTeX-style label
        end
        
        function [gp, upstr] = processInputs(gp, model, varargin)
            % Process optional input parameters
            opt = struct('Db', [], 'upwind_name', 'PhaseUpwindFlag');
            opt = merge_options(opt, varargin{:});
            
            % Set diffusion coefficient if provided
            if ~isempty(opt.Db)
                gp.Db = opt.Db;
            end
            
            % Set upwind name
            gp.upwind_name = opt.upwind_name;
            
            % Create upwind discretization if not provided
            if nargin < 3
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
        end

        function fluxbact = evaluateOnDomain(prop, model, state)
            % Evaluate bacterial diffusive flux
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   fluxbact - Bacterial diffusive flux [cells/(m^2·s)]
            
            % Get upwind flags
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            
            % Get bacterial concentration and max capacity
            nbact = model.getProp(state, 'nbact');
            n0 = model.nbactMax;
                       
            % Get phase properties
            rho = prop.getEvaluatedDependencies(state, 'Density');
            s = model.getProp(state, 's');                
            L_ix = model.getLiquidIndex();

            % Extract liquid phase properties
            if iscell(s)
                sL = s{L_ix};
                rhoL = rho{L_ix};
            else
                sL = s(:, L_ix);
                rhoL = rho(:, L_ix);
            end
            
            % Calculate effective volume fraction
            if iscell(sL)
                Voln = sL{1} .* rhoL{1};
            else
                Voln = sL .* rhoL;
            end
            
            % Compute gradient of bacterial concentration
            gradn = model.operators.Grad(nbact);
            
            % Calculate effective diffusion coefficient
            D_eff = prop.Db .* Voln;
            
            % Compute diffusive flux using upwind scheme
            fluxbact = -prop.faceUpstream(model, state, flag{L_ix}, D_eff) .* gradn;
        end
    end
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Science Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}