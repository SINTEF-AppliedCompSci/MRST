classdef GenericOverallCompositionModel < OverallCompositionCompositionalModel & ExtendedReservoirModel
    properties
        
    end
    
    methods
        function model = GenericOverallCompositionModel(varargin)
            model = model@OverallCompositionCompositionalModel(varargin{:});
            model.OutputProperties = {'ComponentTotalMass'};
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@ReservoirModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Discretize
            [eqs, flux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Define helper variables
            wat = model.water;
            ncomp = numel(names);
            n_hc = ncomp - wat;
            % Assemble equations and add in sources
            
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
            comps = cellfun(@(x, y) {x, y}, X(:, 1+model.water), X(:, 2+model.water), 'UniformOutput', false);
            
            
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                             pressures, sat, mob, rho, ...
                                                             {}, comps, ...
                                                             drivingForces);
            
            for i = 1:numel(eqs)
                if ~isempty(src.cells)
                    eqs{i}(src.cells) = eqs{i}(src.cells) - src.value{i};
                end
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            if false
                cMass = value(model.getProp(state0, 'ComponentTotalMass')');
                massT = sum(cMass, 2);
                scale = (dt./massT);
                for i = 1:ncomp
                    eqs{i} = eqs{i}.*scale;
                end
            else
                massT = model.getComponentScaling(state0);
                scale = (dt./model.operators.pv)./massT;
                for i = 1:ncomp
                    eqs{i} = eqs{i}.*scale;
                end
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Finally assemble
            eqs = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];
        end
        
        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end

        function [state, report] = updateState(model, state, problem, dx, forces)
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end
        
        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'ExtendedFacilityModel')
                model.FacilityModel = ExtendedFacilityModel(model);
            end
            if isempty(model.Components)
                f = model.EOSModel.fluid;
                names_hc = f.names;
                n_hc = numel(names_hc);
                if model.water
                    names = ['water', names_hc];
                else
                    names = names_hc;
                end
                nc = numel(names);
                model.Components = cell(1, nc);
                p = model.FacilityModel.pressure;
                T = model.FacilityModel.T;
                for ci = 1:nc
                    name = names{ci};
                    switch name
                        case 'water'
                            c = ImmiscibleComponent('water', 1);
                        otherwise
                            z = zeros(1, n_hc);
                            z(strcmp(names_hc, name)) = 1;
                            L = standaloneFlash(p, T, z, model.EOSModel);
                            if model.water
                                frac = [0, L, 1-L];
                            else
                                frac = [L, 1-L];
                            end
                            c = EquationOfStateComponent(names{ci}, ci, frac);
                    end
                    model.Components{ci} = c;
                end
            end
            if isempty(model.FlowPropertyFunctions)
                model.FlowPropertyFunctions = CompositionalFlowPropertyFunctions(model);
            end
            model = validateModel@OverallCompositionCompositionalModel(model, varargin{:});
            model.FluxDiscretization.GravityPotentialDifference.saturationWeighting = true;
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.outputFluxes
                f = model.getProp(state, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = [f{:}];
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state, before a possible
            % initialization as AD.'
            [p, sW, z] = model.getProps(state, 'pressure', 'water', 'z');
            z_tol = model.EOSModel.minimumComposition;
            z = ensureMinimumFraction(z, z_tol);
            z = expandMatrixToCell(z);
            cnames = model.EOSModel.fluid.names;
            names = [{'pressure'}, cnames(1:end-1)];
            vars = [p, z(1:end-1)];
            if model.water
                names = [names, {'water'}];
                vars = [vars, {sW}];
            end
            origin = cell(1, numel(names));
            [origin{:}] = deal(class(model));
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state.wellSol);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

        function state = initStateAD(model, state, vars, names, origin)
            isP = strcmp(names, 'pressure');
            state = model.setProp(state, 'pressure', vars{isP});
            removed = isP;
            
            cnames = model.EOSModel.fluid.names;
            ncomp = numel(cnames);
            z = cell(1, ncomp);
            z_end = 1;
            for i = 1:ncomp
                name = cnames{i};
                sub = strcmp(names, name);
                if any(sub)
                    z{i} = vars{sub};
                    z_end = z_end - z{i};
                    removed(sub) = true;
                else
                    fill = i;
                end
            end
            z{fill} = z_end;
            state = model.setProp(state, 'components', z);
            if isa(state.pressure, 'ADI') || isa(z{1}, 'ADI')
                [state.x, state.y, state.L] = model.EOSModel.getPhaseFractionAsADI(state, state.pressure, state.T, state.components);
            end
            if ~isempty(model.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            if model.water
                isWater = strcmp(names, 'water');
                sW = vars{isWater};
                removed(isWater) = true;
                offset = 1;
            else
                offset = 0;
            end
            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
            
            % Now that props have been set up, we can get the saturations
            Z = model.getProps(state, 'PhaseCompressibilityFactors');
            Z_L = Z{offset+1};
            Z_V = Z{offset+2};
            
            L = state.L;
            volL = L.*Z_L;
            volV = (1-L).*Z_V;
            volT = volL + volV;
            sL = volL./volT;
            sV = volV./volT;
            
            if model.water
                void = 1 - sW;
                sL = sL.*void;
                sV = sV.*void;
                s = {sW, sL, sV};
            else
                s = {sL, sV};
            end
            state = model.setProp(state, 's', s);
        end
        
        
        function [sO, sG] = setMinimumTwoPhaseSaturations(model, state, sW, sO, sG, pureVapor, pureLiquid)
            stol = 1e-8;
            if model.water
                sT = sum(state.s, 2);
                if any(pureVapor)
                    sG(pureVapor) = sT(pureVapor) - sW(pureVapor);
                    if isa(sG, 'ADI')
                        sG.val(pureVapor) = max(sG.val(pureVapor), stol);
                    else
                        sG(pureVapor) = max(sG(pureVapor), stol);
                    end
                end

                if any(pureLiquid)
                    sO(pureLiquid) = sT(pureLiquid) - sW(pureLiquid);
                    if isa(sO, 'ADI')
                        sO.val(pureLiquid) = max(sO.val(pureLiquid), stol);
                    else
                        sO(pureLiquid) = max(sO(pureLiquid), stol);
                    end
                end
            end
        end
        
        function forces = validateDrivingForces(model, forces)
            forces = validateDrivingForces@OverallCompositionCompositionalModel(model, forces);
            if ~isempty(forces.W) && numel(forces.W) > 0
                assert(~isempty(model.FacilityModel), ...
                    'FacilityModel must be set up before validating driving forces for a compositional problem with wells!');
                T = model.FacilityModel.T;
                p = model.FacilityModel.pressure;
                wat = model.water;
                assert(isfield(forces.W, 'components'), ...
                    'Wells must have field .components for a compositional model.');
                eos = model.EOSModel;
                for i = 1:numel(forces.W)
                    z = forces.W(i).components;
                    [L, x, y, Z_L, Z_V, rhoL, rhoV] = standaloneFlash(p, T, z, eos);
                    if false
                        % Use flash
                        [sL, sV] = eos.computeSaturations(rhoL, rhoV, x, y, L, Z_L, Z_V);
                        % compi is a mass-fraction in practice
                        L_mass = sL.*rhoL./(sL.*rhoL + sV.*rhoV);
                        comp = [L_mass, 1-L_mass];
                    else
                        % Use the pre-computed definition of light/heavy
                        % components to determine "compi"
                        Z = eos.getMassFraction(z);
                        isEOS = cellfun(@(x) isa(x, 'EquationOfStateComponent'), model.Components);
                        val = cellfun(@(x) x.surfacePhaseMassFractions, model.Components(isEOS), 'UniformOutput', false)';
                        val = vertcat(val{:});
                        val = val(:, (1+wat):end);
                        comp = sum(bsxfun(@times, val, Z'), 1);
                    end
                    if wat
                        assert(~isempty(forces.W(i).compi), ...
                            'W.compi must be present for compositional flow with water phase.');
                        sW = forces.W(i).compi(1);
                        comp = [sW, comp.*(1-sW)];
                        rho = [model.fluid.rhoWS, rhoL, rhoV];
                    else
                        rho = [rhoL, rhoV];
                    end
                    forces.W(i).compi = comp;
                    forces.W(i).rhoS = rho;
                end
            end
        end
    end
end