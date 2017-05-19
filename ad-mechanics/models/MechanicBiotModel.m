classdef MechanicBiotModel < PhysicalModel
    % Two phase oil/water system without dissolution
    properties
        mech;
        rock;
        alpha_scaling;
        S;
        ilu_tol;
    end

    methods
        function model = MechanicBiotModel(G, rock, mech_problem, varargin)
            opt = struct('InputModel', []);
            [opt, rest] = merge_options(opt, varargin{:});

            model = model@PhysicalModel(G, 'stepFunctionIsLinear', true, rest{:});


            % Saving enhanced grid structure within model
            if isempty(opt.InputModel)
                model.G = mrstGridWithFullMappings(model.G);
                model.G = computeGeometryCalc(model.G);
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
                alpha_scaling = 1;
                S             = [];
                ilu_tol       = 1e-4;
                operators = setupOperatorsVEM(model.G, ...
                                              model.mech.C, ...
                                              model.mech.el_bc, ...
                                              model.mech.load, ...
                                              alpha_scaling, S, ilu_tol);

            else
                operators = opt.InputModel.operators;
            end
            model.operators.mech  = operators.mech;
            model.operators.extra = operators.extra;

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose'       , mrstVerbose , ...
                         'reverseMode'   , false       , ...
                         'scaling'       , []          , ...
                         'resOnly'       , false       , ...
                         'history'       , []          , ...
                         'iteration'     , -1          , ...
                         'stepOptions'   , []          , ...
                         'addflux'       , false); % Compatibility only

            opt = merge_options(opt, varargin{:});

            fluidp = drivingForces.fluidp;

            % To be fixed!
            % fbc = drivingForces.fbc;
            fbc = 0;

            xd = model.getProps(state, 'xd');

            if ~opt.resOnly,
                % ADI variables needed since we are not only computing residuals.
                xd = initVariablesADI(xd);
            end

            eqs = equationsMechanicBiot(xd, fluidp, model.G, model.rock, model.operators);

            primaryVars = {'xd'};
            names = {'disp'};
            types = {'disp_dofs'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            problem.iterationNo = opt.iteration;

        end

        function forces = getValidDrivingForces(model)
            % fluid pressure in the volume
            forces.fluidp = [];
        end

        function [fn, index] = getVariableField(model, name)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
            switch(lower(name))
              case {'xd'}
                fn = 'xd';
                index = 1;
              case {'uu'}
                fn = 'uu';
                index = ':';
              case {'u'}
                fn = 'u';
                index = ':';
              case {'stress'}
                fn = 'stress';
                index = ':';
              case {'strain'}
                fn = 'strain';
                index = ':';
              case {'vdiv'}
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
            % add extra model states things from mechanics
            state = addDerivedQuantities(model, state);
        end


        function model = setupOperators(model, alpha_scaling, S, ilu_tol)

            [~, extra] = VEM_linElast(model.G                    , ...
                                      model.mech.C                   , ...
                                      model.mech.el_bc               , ...
                                      model.mech.load                , ...
                                      'alpha_scaling', alpha_scaling , ...
                                      'S', S                         , ...
                                      'linsolve', @(A, rhs) 0 * rhs);

            operators = extra.disc;

            vdiv    = VEM_div(model.G);
            [C, invC, invCi]       = Enu2C(model.mech.Ev, model.mech.nuv, model.G);
            [~, op] = VEM_mrst_vec(model.G, C);%, 'blocksize', model.GMech.cells.num/10);
            strain  = op.WC' * op.assemb';
            stress  = op.D * strain; % The stress will be using "pseudo Voigt's"
                                     % notation, in the sense that it gets an extra
                                     % coefficient 2 for the off-diagonal terms.
            
            % computing and storing useful preconditioner (incomplete Choleski)
            if isnan(ilu_tol)
                fprintf('Skipping ilu\n');
                iL = [];
            else
                iL = shiftedIChol(model.operators.mech.A, 0.5, 'droptol', ilu_tol, 'type', 'ict');
            end

            model.operators.extra = struct('vdiv', vdiv, 'stress', stress, 'strain', ...
                                           strain, 'precond', iL);

        end

    end
end
