classdef MechWaterModel < MechFluidModel


    methods
        function model = MechWaterModel(G, rock, fluid, mech_problem, varargin)

            model = model@MechFluidModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});

        end

        function fluidModel = setupFluidModel(model)
            fluidModel = WaterFluidModel(model.G, model.rock, ...
                                               model.fluid);
        end


        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, wellSol, xd] = model.getProps(state, 'pressure', 'wellSol', ...
                                                     'xd');
            % Properties at previous timestep
            [p0, xd0] = model.getProps(state0, 'pressure', 'xd');

            %Initialization of primary variables

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel; % shortcuts


            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            [p, wellVars{:}, xd] = initVariablesADI(p, wellVars{:}, xd);

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd);

            [w_eqs, w_eqsnames, w_eqstypes, state] = equationsWaterMech(state0, ...
                                                              p, wellVars, ...
                                                              state, fluidModel, ...
                                                              dt, mechTerm, ...
                                                              drivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(w_eqs, mech_eqs);
            names = {w_eqsnames{:}, mech_eqsnames{:}};
            types = {w_eqstypes{:}, mech_eqstypes{:}};

            primaryVars = {'pressure', wellVarNames{:}, 'xd'};
            % make sure that the primary variables defined here match with
            % those of WaterFluidModel and MechanicalModel.

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
