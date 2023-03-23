classdef DualContinuumMechanicModel < MechanicModel
%
% Model that can handle a dual-continuum mechanical system. The inputs can be standard
% boundary conditions and body forces for mechanics, but also a fluid
% pressure, which enters the mechanical equation as grad(pressure).
%
% Even if the equations are linear, we will use the standard nonlinear solver
% to solve the system. It does not make any significant difference as the 
% solver will converge in one Newton step.
%
% SYNOPSIS:
%   model = DualContMechanicModel(G, rock, rock_matrix, mech_problem, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure for the fracture phase
%   rock_matrix  - Rock structure for the matrix
%   mech_problem - Struct containing parameters required to instatiate
%                  mechanical problem
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: DualContMechWaterModel, DualContMechFluidFixedStressSplitModel.m
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
        rock_matrix;  
        constitutive_coefficients_object;
    end

    methods
        function model = DualContinuumMechanicModel(G, rock, rock_matrix, mech_problem, varargin)
        % Constructor for DualContinuumMechanicModel
            model = model@MechanicModel(G, rock, mech_problem, varargin{:});
            
            % Physical properties of matrix rock
            model.rock_matrix  = rock_matrix;
            
            % Create required matrix and fracture stiffness tenors
            coefficient_model_type = evalCoefficientModelType(model);          
            switch coefficient_model_type
                case 'requireInformation'
                    disp('Missing matrix mechanical properties')
                case 'anisotropicCoefficientModels'
                    coefficient_model_handle = str2func(coefficient_model_type);
                    model.constitutive_coefficients_object = coefficient_model_handle(model);
                    [model.mech.C, model.mech.invC, model.mech.invCi] = ...
                            AnisoEnu2C(model.mech.E, model.mech.nu, model.mech.mu, model.G);
                case 'isotropicStiffCoefficientModels'
                    coefficient_model_handle = str2func(coefficient_model_type);
                    model.constitutive_coefficients_object = coefficient_model_handle(model);
                case 'isotropicVoidCoefficientModels'
                    disp(['No stiffness properties given for the fracture phase, ',...
                         'using void space coefficient models']);
                    coefficient_model_handle = str2func(coefficient_model_type);
                    model.constitutive_coefficients_object = coefficient_model_handle(model);
            end
                
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            % Assemble the equations for the mechanical system
            opt = struct('Verbose'       , mrstVerbose , ...
                         'reverseMode'   , false       , ...
                         'scaling'       , []          , ...
                         'resOnly'       , false       , ...
                         'history'       , []          , ...
                         'iteration'     , -1          , ...
                         'stepOptions'   , []          , ...
                         'addflux'       , false); % Compatibility only

            opt = merge_options(opt, varargin{:});
            
            % Shortcut for constitutive model
            constitutiveModel = model.constitutive_coefficients_object;
                        
            % The fluid pressure stimulates the mechanical system. It is
            % given as a driving force.
            fluidp = drivingForces.fluidp;
            fluidp_matrix = drivingForces.fluidp_matrix;
            
            xd = model.getProps(state, 'xd');

            if ~opt.resOnly
                xd = initVariablesADI(xd);
            end

            eqs = equationsDCPoroMechanics(xd, model, fluidp, fluidp_matrix);

            primaryVars = {'xd'};
            names = {'disp'};
            types = {'disp_dofs'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            problem.iterationNo = opt.iteration;
        end

        function forces = getValidDrivingForces(model)
        % The fluid pressure stimulates the mechanical system. It is 
        % given as a driving force.
            forces.fluidp = [];
            forces.fluidp_matrix = [];
        end

        function [fn, index] = getVariableField(model, name, varargin)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
            switch(lower(name))
              case {'uu'}
                % Displacement field as a matrix (one column per Cartesian direction)
                fn = 'uu';
                index = ':';
              case {'u'}
                % Displacement field given as a column vector where the cartesian components are
                % stabbed.
                fn = 'u';
                index = ':';
              case {'xd'}
                % Same as 'u' but the degree of freedom where the Dirichlet conditions (fixed
                % displacement) are removed
                fn = 'xd';
                index = 1;
%               case {'local_strain'}
%                 % Local volumetric strain 
%                 fn = 'local_strain';
%                 index = ':';
              case {'stress'}
                % Stress tensor (Voigt notation, one column per component)
                fn = 'stress';
                index = ':';
              case {'strain'}
                % Strain tensor (Voigt notation, one column per component)
                fn = 'strain';
                index = ':';
              case {'vol_strain_frac'}
                % Intrinsic volumetrix fracture strain
                fn = 'vol_strain_frac';
                index = ':';
              case {'vol_strain_mat'}
                % Intrinsic volumetrix matrix strain
                fn = 'vol_strain_mat';
                index = ':';
              case {'vdiv'}
                % Volume weighted divergence field of displacement.
                fn = 'vdiv';
                index = ':';
              otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@MechanicModel(model, name);
            end
        end

        function [primaryVars, fds] = getAllVarsNames(model)
        % List the variables that are used by the model. Used when 
        % coupling the mechanical model with a fluid model.
            primaryVars = {'xd'};
            fds = {'xd', 'uu', 'u', 'local_strain', 'stress', 'strain', ...
                   'strain_frac', 'strain_mat', 'vdiv'};
        end
        
        function [coefficient_model_type] = evalCoefficientModelType(model)
        % Function to evaluate the type of coefficient model that we need
        % to use. Returns a string, which will be converted to a function
        % handle used to identify the appropriate models to be used in the
        % calculation of the coefficient models.
            if ~isfield(model.mech, 'E_m') && ~isfield(model.mech, 'nu_m')
                coefficient_model_type = 'requireInformation';     
            elseif isfield(model.mech, 'E_f') && isfield(model.mech, 'nu_f')
                try
                % test for anisotropy by ensuring we have the required
                % number of mechanical parameters to construct the
                % stiffness and compliance tensors
                    assert((size(model.mech.E_f, 2)) == model.G.griddim)
                    coefficient_model_type = 'anisotropicCoefficientModels';
                catch
                    coefficient_model_type = 'isotropicStiffCoefficientModels';
                end
            else                  
                coefficient_model_type = 'isotropicVoidCoefficientModels'; 
            end
        end
        
    end
end
