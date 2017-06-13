classdef MechanicModel < PhysicalModel
% Model that can handle a mechanical system. The inputs can be standard
% boundary conditions and body forces for mechanics, but also a fluid
% pressure, which enters the mechanical equation as grad(pressure).
%
% Even if the equations are linear, we will use a standard nonlinear solver
% (which should converge in one Newton step).
%
% SYNOPSIS:
%   model = MechanicModel(G, mech_problem, varargin)
%
% DESCRIPTION:
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
%
% SEE ALSO:
%

    
    % Mechanical model
    properties
        % Structure that contains the mechanical parameter of the system
        mech; 
        % Structure that contains rock properties (Biot coefficient)
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
                    Enu2C(model.mech.Ev, model.mech.nuv, model.G);
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
            forces.fluidp = [];
        end

        function [fn, index] = getVariableField(model, name)
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
                index = 1;
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
                [fn, index] = getVariableField@PhysicalModel(model, name);
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
