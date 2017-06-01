classdef BlackOilFixedStressFluidModel < ThreePhaseBlackOilModel

    properties
        pressCoef;
    end

    methods
        function model = BlackOilFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            model.disgas = true;
            model.vapoil = false;
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only

            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, sW, sG, rs, rv, wellSol] = model.getProps(state, 'pressure', ...
                                                                 'water', ...
                                                                 'gas', ...
                                                                 'rs', 'rv', ...
                                                                 'wellSol');
            % Properties at previous timestep
            [p0, sW0, sG0, rs0, rv0] = model.getProps(state0, 'pressure', ...
                                                              'water', ...
                                                              'gas', 'rs', 'rv');

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
                [p, sW, x, wellVars{:}] = ...
                    initVariablesADI(p, sW, x, wellVars{:});
            end

            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;
            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;

            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');

            [eqs, names, types, state] = equationsBlackOilMech(state0, st0, ...
                                                              p, sW, x,  rs, ...
                                                              rv, st, wellVars, ...
                                                              state, model, ...
                                                              dt, mechTerm, ...
                                                              otherDrivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            primaryVars = {'pressure', 'sW', gvar, 'qWs', 'qOs', 'qGs', 'bhp'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ThreePhaseBlackOilModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end

        function fds = getAllVarsNames(model)
            fds = {'wellSol', 'pressure', 's', 'rs', 'rv'};
        end


    end

end
