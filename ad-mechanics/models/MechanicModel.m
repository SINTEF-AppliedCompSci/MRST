classdef MechanicModel < PhysicalModel
%
% Model that can handle a mechanical system. The inputs can be standard
% boundary conditions and body forces for mechanics, but also a fluid
% pressure, which enters the mechanical equation as grad(pressure).
%
% Even if the equations are linear, we will use the standard nonlinear
% solver to solve the system. It does not make any significant
% difference as the the solver will converge in one Newton step.
%
% SYNOPSIS:
%   model = MechanicModel(G, rock, mech_problem, varargin)
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechBlackOilModel, MechOilWaterModel, MechWaterModel,
%   MechFluidFixedStressSplitModel

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
        % Structure that contains the mechanical parameter of the system
        % This structure should contain the fields:
        %
        % E     - Young's modulus (one entry per cell)
        % nu    - Poisson ratio (one entry per cell)
        % el_bc - Structure describing the boundary condition (see VEM_linElast)
        % load  - Structure giving the volumetric forces (see VEM_linElast)
        %
        mech;
        % Structure that contains the rock properties. The structure must
        % contain the Biot coefficient
        %
        % alpha - Biot coefficient (one entry per cell)
        %
        rock;

        % Parameters that are used in the VEM assembly, see setupOperatorsVEM and
        % VEM_linElast.
        alpha_scaling;
        S;
        ilu_tol;

    end

    methods
        function model = MechanicModel(G, rock, mech_problem, varargin)
        % Constructor for MechanicModel

            opt = struct('InputModel', []);
            [opt, rest] = merge_options(opt, varargin{:});

            model = model@PhysicalModel(G, rest{:});

            % Process the grid for mechanical computation
            if any(strcmpi('createAugmentedGrid', model.G.type))
                model.G = createAugmentedGrid(model.G);
            end

            % Physical properties of rock and fluid
            model.mech  = mech_problem;
            model.rock  = rock;

            % Compute stiffness tensor C, if not given
            if ~isfield(model.mech, 'C')
                [model.mech.C, model.mech.invC, model.mech.invCi] = ...
                    Enu2C(model.mech.E, model.mech.nu, model.G);
            end

            if isempty(opt.InputModel)
                alpha_scaling = 1;  % default values
                S             = []; % default values
                operators = setupOperatorsVEM(model.G, ...
                                              model.mech.C, ...
                                              model.mech.el_bc, ...
                                              model.mech.load, ...
                                              alpha_scaling, S);

            else
                operators = opt.InputModel.operators;
            end
            model.operators = operators;

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

            % The fluid pressure stimulates the mechanical system. It is
            % given as a driving force.
            fluidp = drivingForces.fluidp;


            xd = model.getProps(state, 'xd');

            if ~opt.resOnly,
                xd = initVariablesADI(xd);
            end

            eqs = equationsPoroMechanics(xd, model, fluidp);

            primaryVars = {'xd'};
            names = {'disp'};
            types = {'disp_dofs'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            problem.iterationNo = opt.iteration;

        end

        function forces = getValidDrivingForces(model)
            % The fluid pressure stimulates the mechanical system. It is given as a driving
            % force.
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.fluidp = [];
        end

        function [fn, index] = getVariableField(model, name, varargin)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
            switch(lower(name))
              case {'uu'}
                % displacement field as a matrix (one column per Cartesian direction)
                fn = 'uu';
                index = ':';
              case {'u'}
                % displacement field given as a column vector where the cartesian components are
                % stabbed.
                fn = 'u';
                index = ':';
              case {'xd'}
                % same as 'u' but the degree of freedom where the Dirichlet conditions (fixed
                % displacement) are removed
                fn = 'xd';
                index = ':';
              case {'stress'}
                % Stress tensor (Voigt notation, one column per component)
                fn = 'stress';
                index = ':';
              case {'strain'}
                % Strain tensor (Voigt notation, one column per component)
                fn = 'strain';
                index = ':';
              case {'vdiv'}
                % volume weighted divergence field of displacement.
                fn = 'vdiv';
                index = ':';
              otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@PhysicalModel(model, state, problem, dx, drivingForces);
            % add extra model states  from mechanics
            state = addDerivedQuantities(model, state);
        end

        function [primaryVars, fds] = getAllVarsNames(model)
        % list the variables that are used by the model. Used when coupling
        % the mechanical model with a fluid model.
            primaryVars = {'xd'};
            fds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
        end

    end
end
