classdef GenericBlackOil < ThreePhaseBlackOilModel & ExtendedReservoirModel
    properties
        
    end
    
    methods
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            opt = struct('Verbose',     mrstVerbose,...
                        'reverseMode', false,...
                        'resOnly',     false,...
                        'iteration',   -1, ...
                        'drivingForces0', []);
            opt = merge_options(opt, varargin{:});


            % Define primary variables
            if opt.reverseMode
                [state0, primaryVars] = model.getReverseStateAD(state0);
                % The model must be validated with drivingForces so that the
                % FacilityModel gets updated.
                model = model.validateModel(drivingForces);
                state = model.getStateAD(state, false);
            else
                [state, primaryVars] = model.getStateAD(state, ~opt.resOnly);
            end


            [eqs, divTerms, names, types] = conservationEquationsBlackOil(state0, state, model, dt, drivingForces);

            dissolved = model.getDissolutionMatrix(state.rs, state.rv);
            % Add in and setup well equations

            wellVars = state.FacilityState.primaryVariables;
            w = state.FacilityState.names;
            nw = numel(state.wellSol);
            wellMap = struct('isBHP', strcmp(w, 'bhp'), 'isRate', ...
                              ismember(w, {'qOs', 'qGs', 'qWs'}), ...
                             'extraMap', zeros(nw, nw));

            p = state.pressure;
            mob = state.FlowProps.Mobility;
            rho = state.FlowProps.Density;
            sat = state.s;
            pressures = state.FlowProps.PhasePressures;


            [eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                             pressures, sat, mob, rho, ...
                                                             dissolved, {}, ...
                                                             drivingForces);
            [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                              names, types, ...
                                                              state0.wellSol, ...
                                                              state.wellSol, ...
                                                              wellVars, wellMap, ...
                                                              p, mob, rho, dissolved, ...
                                                              {}, dt, opt);
            for i = 1:numel(divTerms)
                eqs{i} = eqs{i} + divTerms{i};
            end
            state = value(state);
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end

        
        
        function c = getComponentDensity(model, state, name, ph)
            if nargin < 4
                ph = model.getPhaseNames();
            end
            nph = numel(ph);

            c = cell(nph, 1);
            switch name
                case 'oil'
                    oix = strcmpi(ph, 'o');
                    if model.vapoil
                        [b, rv] = model.getProps(state, 'ShrinkageFactors', 'rv');
                        rhoS = model.getSurfaceDensities();
                        gix = strcmpi(ph, 'g');
                        rhoOS = rhoS(oix);
                        bO = b{oix};
                        bG = b{gix};
                        c{oix} = rhoOS.*bO;
                        c{gix} = rv.*rhoOS.*bG;
                    else
                        c{oix} = 1;
                    end
                case 'gas'
                    error();
                otherwise
                    c = getComponentDensity@ExtendedReservoirModel(model, state, name, ph);
            end
        end
    end
end