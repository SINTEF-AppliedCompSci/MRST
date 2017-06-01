classdef MechBlackOilModel2 < MechFluidModel

    properties
    end

    methods
        function model = MechBlackOilModel2(G, rock, fluid, mech_problem, varargin)

            model = model@MechFluidModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});

        end

        function fluidModel = setupFluidModel(model)
            fluidModel = BlackOilFluidModel(model.G, model.rock, ...
                                                 model.fluid);
            fluidModel.disgas = true;
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

            %Initialization of primary variables

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel; % shortcuts

            st  = fluidModel.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
            st0 = fluidModel.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
            if fluidModel.disgas || fluidModel.vapoil
                % X is either Rs, Rv or Sg, depending on each cell's saturation status
                x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
                gvar = 'x';
            else
                x = sG;
                gvar = 'sG';
            end

            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, x, wellVars{:}, xd] = initVariablesADI(p, sW, x, ...
                                                               wellVars{:}, ...
                                                               xd);
            end

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd);

            [bo_eqs, bo_eqsnames, bo_eqstypes, state] = equationsBlackOilMech(state0, ...
                                                              st0, p, sW, x, ...
                                                              rs, rv, st, ...
                                                              wellVars, state, ...
                                                              fluidModel, dt, ...
                                                              mechTerm, ...
                                                              drivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics2(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(bo_eqs, mech_eqs);
            names = {bo_eqsnames{:}, mech_eqsnames{:}};
            types = {bo_eqstypes{:}, mech_eqstypes{:}};

            primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}, 'xd'};

            model.fluidModel.primaryVars = {'pressure', 'sW', gvar, ...
                                wellVarNames{:}};
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                              xd0, p, xd)

            opmech = model.mechModel.operators.mech;
            fluidp = p;
            mechTerm.new = opmech.div*xd;
            mechTerm.old = opmech.div*xd0;

        end


    end
end
