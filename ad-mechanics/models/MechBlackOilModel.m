classdef MechBlackOilModel < ThreePhaseBlackOilModel
    properties
        mech;
        alpha_scaling;
        S;
    end

    methods
        function model = MechBlackOilModel(G, rock, fluid, mech_problem, varargin)

            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model.disgas = true;
            
            model.alpha_scaling = 1;
            model.S = [];
            model = merge_options(model, varargin{:});

            model.G = createAugmentedGrid(model.G);
            model.mech = mech_problem;

            % Compute stiffness tensor C, if not given
            if ~isfield(model.mech, 'C')
                model.mech.C = Enu2C(model.mech.Ev, model.mech.nuv, model.G);
            end
            % Blackoil -> use CNV style convergence
            model.useCNVConvergence = false;

            operators = setupOperatorsVEM(model.G, ...
                                          model.mech.C, ...
                                          model.mech.el_bc, ...
                                          model.mech.load, ...
                                          model.alpha_scaling, ...
                                          model.S);

            model.operators.mech  = operators.mech;
            model.operators.extra = operators.extra;

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, sW, sG, rs, rv, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                                     'water', ...
                                                                     'gas', ...
                                                                     'rs', ...
                                                                     'rv', ...
                                                                     'wellSol', ...
                                                                     'xd');
            % Properties at previous timestep
            [p0, sW0, sG0, rs0, rv0, xd0] = model.getProps(state0, 'pressure', ...
                                                                   'water', ...
                                                                   'gas', 'rs', ...
                                                                   'rv', ...
                                                                   'xd');
            
            %Initialization of primary variables ----------------------------------
            st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
            st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
            if model.disgas || model.vapoil
                % X is either Rs, Rv or Sg, depending on each cell's saturation status
                x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
                gvar = 'x';
            else
                x = sG;
                gvar = 'sG';
            end

            [wellVars, wellVarNames, wellMap] = ...
                model.FacilityModel.getAllPrimaryVariables(wellSol);
            
            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, x, wellVars{:}, xd] = initVariablesADI(p, ...
                                                                  sW, x, ...
                                                                  wellVars{:}, ...
                                                                  xd);
            end

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd);

            [blackoileqs, state] = equationsBlackOilMech(state0, st0, p, sW, ...
                                                         x, rs, rv, st, ...
                                                         wellVars, state, ...
                                                         model, dt, mechTerm, ...
                                                         drivingForces, ...
                                                         'iteration', ...
                                                         opt.iteration);
            mecheqs = equationsPoroMechanics(xd, fluidp, model.G, model.rock, ...
                                            model.operators);

            eqs = horzcat(blackoileqs, mecheqs);
            primaryVars = {'pressure', 'sW', gvar, 'qWs', 'qOs', 'qGs', 'bhp', ...
                          'xd'};
            names = {'water', 'oil', 'gas', wellVarNames{:}, 'disp'};
            types = {'cell', 'cell', 'cell', 'perf', 'perf', 'perf', 'well', ...
                     'disp_dofs'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


        function [fn, index] = getVariableField(model, name)
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
                [fn, index] = getVariableField@ThreePhaseBlackOilModel(model, name);
            end
        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                              xd0, p, xd)
            
            s = model.operators;
            fluidp = p;
            mechTerm.new = s.mech.div*xd;
            mechTerm.old = s.mech.div*xd0;

        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, problem, dx, drivingForces);
            % add extra model states things from mechanics
            state = addDerivedQuantities(model,state);
        end

        function model = setupOperators(model, G, rock, varargin)

            
            % Set up divergence / gradient / transmissibility operators for flow
            model = setupOperators@ReservoirModel(model, G, rock, varargin{:});

            operators = setupOperatorsVEM(model.G, model.mech.el_bc, ...
                                                   model.mech.load, ...
                                                   model.alpha_scaling, ...
                                                   model.S);
            model.operators.mech = operators.mech;
            model.operators.extra = operators.extra;
            
        end

    end
end
