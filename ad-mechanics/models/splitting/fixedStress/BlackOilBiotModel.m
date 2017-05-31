classdef BlackOilBiotModel < ThreePhaseBlackOilModel

    properties
        fluidModelType;
    end


    methods
        function model = BlackOilBiotModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            fluidModelType = 'blackoil';
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

            bhp    = vertcat(wellSol.bhp);
            qWs    = vertcat(wellSol.qWs);
            qOs    = vertcat(wellSol.qOs);
            qGs    = vertcat(wellSol.qGs);
            
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
            
            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, x, qWs, qOs, qGs, bhp] = ...
                    initVariablesADI(p, sW, x, qWs, qOs, qGs, bhp);
            end

            mechTerm = drivingForces.divTerm;
            otherDrivingForces = rmfield(drivingForces, 'divTerm');


            [eqs, state] = equationsBlackOilBiot(state0, st0, p, sW, ...
                                                         x, bhp, qWs, qOs, qGs, rs, ...
                                                         rv, st, state, model, ...
                                                         dt, mechTerm, ...
                                                         drivingForces, ...
                                                         'iteration', ...
                                                         opt.iteration);

            primaryVars = {'pressure', 'sW', gvar, 'qWs', 'qOs', 'qGs', 'bhp'};
            names = {'water', 'oil', 'gas' 'waterWells', 'oilWells', 'gasWells', ...
                     'closureWells'};
            types = {'cell', 'cell', 'cell', 'perf', 'perf', 'perf', 'well'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ThreePhaseBlackOilModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.divTerm = [];
        end

    end

end
