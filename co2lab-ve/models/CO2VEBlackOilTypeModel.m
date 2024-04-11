classdef CO2VEBlackOilTypeModel < ReservoirModel & GenericReservoirModel
    
    properties
        hysteresis
        disgas  % the name of this field should not be changed, in order to
                % play nice with the rest of AD-OO.
    end

    methods
        function model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
            
            model = model@ReservoirModel(Gt, varargin{:});
            model.rock = rock2D;
            model.fluid = fluid;
            [model.water, model.oil, model.gas] = deal(true, false, true);
            model.gravity = norm(gravity);
            
            % determine if model includes hysteresis and/or dissolution
            model.hysteresis = model.fluid.res_gas > 0;
            model.disgas = isfield(model.fluid, 'dis_max'); 
            
            model = model.setupOperators(Gt, rock2D);
            model.Components = {MiscibleWaterComponent('brine'), ...
                                CO2Component('CO2', model.disgas)};
        end
        
% ------------------------------------------------------------------------        
        
        function model = setupOperators(model, Gt, rock)
        % Compute vertially-integrated transmissibilities 
            rock_tmp      = rock; 
            rock_tmp.perm = rock.perm .* Gt.cells.H; 
            T             = computeTrans(Gt, rock_tmp); 
            cf            = Gt.cells.faces(:, 1); 
            nf            = Gt.faces.num; 
            T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);       
            
            % Computing vertically-integrated pore-volume
            pv = poreVolume(Gt, rock); 
            pv = pv .* Gt.cells.H; 
            
            model.operators = setupOperatorsTPFA(Gt, rock, 'porv', pv, 'trans', T);
            model.operators.T_all = T;
        end

% ------------------------------------------------------------------------
        
function rhoS = getSurfaceDensities(model, regions, phases)
% getSurfaceDensities@ReservoirModel removes (ADI) derivatives from results, so
% we override it here (to allow optimization that involves surface densities).
    rhoWS = model.fluid.rhoWS;
    rhoGS = model.fluid.rhoGS;
    
    if nargin > 1 && ~isempty(regions)
        rhoWS = rhoWS(regions, :);
        rhoGS = rhoGS(regions, :);
    end
    rhoS = {rhoWS, rhoGS};
end

% ------------------------------------------------------------------------

        function [eqs, names ,types, state] = ...
                getModelEquations(model, state0, state, dt, drivingForces)
            
            [acc, flux, names, types] = ...
                model.FlowDiscretization.componentConservationEquations(...
                    model, state, state0, dt);
            
            % add sources from wells
            src = model.FacilityModel.getComponentSources(state);
            acc = model.insertSources(acc, src);
            
            % assemble equations
            eqs = cell(1, numel(acc));
            for i = 1:numel(acc)
                eqs{i} = model.operators.AccDiv(acc{i}, flux{i});
            end
            pv = model.getProp(state, 'PoreVolume');
            b  = model.getProp(state, 'ShrinkageFactors');
            [bW, bG] = deal(b{:});
            
            if model.hasDissolution()

                [co2dis, co2dis0] = deal(model.dissolvedCO2Mass(state), ...
                                         model.dissolvedCO2Mass(state0));
                
                acc_dis = (co2dis - co2dis0) / dt;
                flux_dis = model.getProp(state, 'DissolvedFlux');
                eq_dis = model.operators.AccDiv(acc_dis, flux_dis);
                
                if model.hasRateDrivenDissolution()
                    % Add conservation equation for dissolved gas
                    eta = model.computeDissolutionTransfer(state0, state, dt);
                    eta_mass = eta .* model.fluid.rhoGS .* bW;
                    eqs = [eqs, {eq_dis - eta_mass}];
                    names = [names, {'dissolution'}];
                    types = [types, {'cell'}];
                    
                    % amount of saturation "eaten up" by dissolution.  Will  be
                    % used to compute depletion of residual staturation further below.
                    s_depletion = eta_mass;
                else
                    assert(model.hasInstantDissolution())
                    
                    s_depletion = max(eq_dis, 0);
                end
            else
                % no dissolution, so no depletion of residual saturation
                s_depletion = 0;
            end
            
            if model.hysteresis
                sG = model.getProp(state, 'sg');
                sGmax = model.getProp(state, 'sGmax');
                sGmax0 = model.getProp(state0, 'sGmax');

                names = [names, {'hysteresis'}];
                types = [types, {'cell'}];

                fac = (1 - model.fluid.res_water) ./ model.fluid.res_gas ./ ...
                      pv ./ model.fluid.rhoGS ./ bG;
                eqs = [eqs, { sGmax - max(sGmax0 - s_depletion .* dt .* fac, sG) }];
                
            end
            
            % get well equations
            if ~isempty(model.FacilityModel)
                [weqs, wnames, wtypes, state] = ...
                    model.FacilityModel.getModelEquations(state0, state, dt, ...
                                                                  drivingForces);
                % concatenate
                eqs = [eqs, weqs];
                names = [names, wnames];
                types = [types, wtypes];
            end
            
            % add in boundary conditions
            [eqs, state, src] = ...
                model.addBoundaryConditionsAndSources(eqs, names, types, ...
                                                      state, drivingForces);
        end

% ------------------------------------------------------------------------        
        function m = dissolvedCO2Mass(model, state)
            w_ph_ix = model.getPhaseIndex('W');
            co2_comp_ix = 2; % @@ hard coded for now
            
            masses = model.getProp(state, 'ComponentPhaseMass');
            m = masses{co2_comp_ix, w_ph_ix};
        end

% ------------------------------------------------------------------------        
        function m = undissolvedCO2Mass(model, state)
            co2_ph_ix = model.getPhaseIndex('G');
            co2_comp_ix = 2; % @@ hard coded for now
            
            masses = model.getProp(state, 'ComponentPhaseMass');
            m = masses{co2_comp_ix, co2_ph_ix};
        end
        
        
% ------------------------------------------------------------------------        

        function eta = computeDissolutionTransfer(model, state0, state, dt)

            % compute maximum possible transfer rate (considering how much is already
            % dissolved, and how much CO2 is available for further dissolution)
            
            [pv, pv0] = deal(model.getProp(state, 'PoreVolume'), ...
                             model.getProp(state0, 'PoreVolume'));
            [b, b0] = deal(model.getProp(state, 'ShrinkageFactors'), ...
                           model.getProp(state0, 'ShrinkageFactors'));
            [sW, sW0] = deal(model.getProp(state, 'sw'), ...
                             model.getProp(state0, 'sw'));
            [sG, sG0] = deal(model.getProp(state, 'sg'), ...
                             model.getProp(state0, 'sg'));
            
            bW = b{1}; bW0 = b0{1};
            bG = b{2}; bG0 = b0{2};
            
            rsmax = model.fluid.dis_max;
            rs0 = model.getProp(state0, 'rs');
            
            state_rsmax = state;
            state_rsmax.rs = state_rsmax.rs * 0 + rsmax;
            state_rsmax = model.initStateFunctionContainers(state_rsmax);
            
            [co2dis_max, co2dis0] = deal(model.dissolvedCO2Mass(state_rsmax), ...
                                         model.dissolvedCO2Mass(state0));
                
            flux_dis = model.getProp(state, 'DissolvedFlux');

            % compute the maximum amount of CO2 that can be absorbed in the
            % current situation over the next timestep.
            max_demand = (co2dis_max - co2dis0) / dt; 
            max_demand = model.operators.AccDiv(max_demand, flux_dis);
            
            % how much dissolved CO2 could maximally be supplied from the
            % remaining CO2 phase in the cell
            max_supply = model.undissolvedCO2Mass(state0) / dt;
            
            src = model.FacilityModel.getComponentSources(state);
            if ~isempty(src.cells)
                max_supply(src.cells) = max_supply(src.cells) + value(src.value{2});
            end
            
            phase_fluxes = model.getProp(state, 'ComponentPhaseFlux');
            co2_phase_flux = phase_fluxes{2, 2}; % @@ second component (CO2) in second phase (CO2)
            
            co2_inflow = -1 * model.operators.C' * co2_phase_flux;
            max_supply = max_supply + max(value(co2_inflow), 0.0);

            % The maximum amount of CO2 mass that can be dissolved in a cell
            % will be limited by the most strict of the two above bounds
            max_transfer = max(min(max_demand, max_supply), 0);

            % convert from mass rate to rs rate
            max_transfer  = max_transfer ./ (bW .* model.fluid.rhoGS);
            
            % The vertical zone where there is both brine and CO2 is considered instantly
            % saturated.  We must ensure that sufficient CO2 is dissolved to
            % fulfil this at all times.
            sW_touched = model.brineSaturationInTwoPhaseZone(state);
            min_dissolved = model.fluid.dis_max .* sW_touched .* bW .* pv; % @@ NB: This assumes 
            cur_dissolved = rs0 .* sW .* bW .* pv;                         % vertical heterogeneity.
            min_transfer = max(min_dissolved - cur_dissolved, 0) / dt;     % Generalize this!
            

            % The amount of actually dissolved CO2 is also limited by the
            % actual dissolution rate
            basic_rate = model.fluid.dis_rate .* pv ./ model.G.cells.H; % rate per area multiplied
                                                                        % by CO2/brine interface
                                                                        % area in cell
            eta =  min(basic_rate, max_transfer);
            %eta = min(max(basic_rate, min_transfer), max_transfer);
            %eta = max(min(basic_rate, max_transfer), min_transfer);
        end

        
% ------------------------------------------------------------------------        
        function sW_touched = brineSaturationInTwoPhaseZone(model, state)
            % @@ This function must be generalized to the vertically non-heterogeneous case!
            sG = model.getProp(state, 'sg');
            sGmax = model.getProp(state, 'sGmax');
            sW = 1 - sG;
            
            [rw, rg] = deal(model.fluid.res_water, model.fluid.res_gas);
            
            [h, h_max] = upscaledSat2height(sG, sGmax, model.G, ...
                                            'resSat', [rw, rg]); % @@ additional arguments for non-heterogeneous
            H = model.G.cells.H;
            
            sW_untouched = (H - h_max) ./ H; % part of water saturation that were always below CO2 plume
            
            sW_touched = sW  - sW_untouched; % part of water saturation that at some point was
                                             % in contact with the CO2
        end
% ------------------------------------------------------------------------        

        function model = validateModel(model, varargin)
            model = validateModel@ReservoirModel(model, varargin{:});
        end

% ------------------------------------------------------------------------        

        function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
        end

% ------------------------------------------------------------------------        

        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@ReservoirModel(model, varargin{:});
        
            % FlowPropertyFunctions
            flowprops = model.FlowPropertyFunctions;
            flowprops = flowprops.setStateFunction('CapillaryPressure', ...
                                                   CO2VECapillaryPressure(model));
            flowprops = flowprops.setStateFunction('RelativePermeability', ...
                                                   CO2VERelativePermeability(model));
            model.FlowPropertyFunctions = flowprops;

            % FluxDiscretization
            flowdisc = model.FlowDiscretization;
            flowdisc = flowdisc.setStateFunction('DissolvedFlux', ...
                                                 CO2VEDissolvedFlux(model));
            model.FlowDiscretization = flowdisc;
            
            % PVTPropertyFunctions
            pvtprops = model.PVTPropertyFunctions;
            pvtprops = pvtprops.setStateFunction('Density', CO2VEBlackOilDensity(model));
            pvtprops = pvtprops.setStateFunction('ShrinkageFactors', ...
                                                 ShrinkageFactors(model,  ...
                                                              'usePhasePressures', false));
            
            model.PVTPropertyFunctions = pvtprops;
        end

% ------------------------------------------------------------------------        
        
        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            
            if model.hasInstantDissolution()

                % In the instant dissolution model, we need to 'unpack' 
                % saturation and 'rs' from the single primary variable 'X'.
                
                % handle variable 'X' separately
                incX = model.getIncrement(dx, problem, 'X');
                rs_orig = model.getProp(state, 'rs');
                unsat = rs_orig < model.fluid.dis_max;

                % update rs field
                state = model.updateStateFromIncrement(state, incX .* unsat, ...
                                                       problem, 'rs', inf, inf);
                state = model.capProperty(state, 'rs', 0, model.fluid.dis_max);
                                                       
                % update saturation
                ds = incX;
                ds(unsat) = 0;
                state = model.updateStateFromIncrement(state, ds, problem, 's', ...
                                                       inf, model.dsMaxAbs);
                state = model.capProperty(state, 's', 0, 1 - model.fluid.res_water);
                
                % ensure that cells with no remaining water saturation 'unsat' by making 
                % all such cells have a rs-value strictly less than fluid.dis_max
                sg_updated = model.getProp(state, 's');
                sg_updated = sg_updated(:, model.getPhaseIndex('G'));
                turned_unsat = (ds < 0) & (sg_updated == 0);
                tiny = zeros(size(sg_updated));
                tiny(turned_unsat) = -1 * eps;
                state = model.updateStateFromIncrement(state, tiny, problem, 'rs', inf, inf);
                
                % Get rid of the 'X' variable
                [problem.primaryVariables, ix] = ...
                    model.stripVars(problem.primaryVariables, {'X'});
                dx(ix) = [];
            end
            
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                         drivingForces);
            
            state.sGmax = max(state.sGmax, state.s(:,2)); % enforce sGmax >= sg
        end


% ------------------------------------------------------------------------        
        
        function [state, report] = ...
                updateAfterConvergence(model, state0, state, dt, drivingForces)
      
            [state, report] = ...
                updateAfterConvergence@ReservoirModel(model, ...
                                                      state0, state, dt, drivingForces);
            
            if model.outputFluxes
                % @@ Can this be done so we won't have to back-up the
                % boundary fluxes?

                % store already-registered boundary fluxes
                if ~isempty(drivingForces.bc)
                    bnd_flux_ixs = drivingForces.bc.face;
                    bnd_fluxes = state.flux(bnd_flux_ixs, :);
                end
                
                % add internal fluxes (@@ will overwrite boundary fluxes)
                qWG = model.getProp(state, 'PhaseFlux');
                state = storeFluxes(model, state, qWG{1}, [], qWG{2});
                
                % copy back the overwritten boundary fluxes
                if ~isempty(drivingForces.bc)
                    state.flux(bnd_flux_ixs, :) = bnd_fluxes;
                end
            end
        end

% ------------------------------------------------------------------------        

        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
              case {'sgmax'}
                index = 1;
                fn = 'sGmax';
              case {'rs'}
                index = 1;
                fn = 'rs';
              otherwise
                [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
            end
        end

% ------------------------------------------------------------------------        
        
       function names = getComponentNames(model)
           names = cellfun(@(c) c.name, model.Components, 'uniformoutput', false);    
           %names = getComponentNames@ReservoirModel(model);
        end

% ------------------------------------------------------------------------        
    
        function [vars, names, origin] = getPrimaryVariables(model, state)
    
            [p, s] = model.getProps(state, 'pressure', 'saturation');
            
            sG = s(:,2); % water saturation.  s(:,2) is CO2 saturation
            
            [vars, names, origin] = deal({p}, {'pressure'}, {class(model)});
            
            if model.hasInstantDissolution()
                
                % Instant dissolution model
                rs = model.getProp(state, 'rs');
                unsat = rs < model.fluid.dis_max;
                % X represents phase saturation in cells where dissolution has
                % reached its maximum, and amount of dissolved CO2 in other cells
                X = sG;
                X(unsat) = rs(unsat);
                vars = [vars, X];
                names = [names, 'X'];
                origin = [origin, class(model)];
                
            elseif model.disgas
                
                % Finite rate dissolution model
                rs = model.getProp(state, 'rs');
                vars = [vars, sG, rs];
                names = [names, 'sG', 'rs'];
                origin = [origin, class(model), class(model)];
                
            else
                
                % We do not consider dissolution, and use sG directly as 
                % primary variable
                vars = [vars, sG];
                names = [names, 'sG'];
                origin = [origin, class(model)];
                
            end
            
            if model.hysteresis
                
                % add sGmax as primary variable
                vars = [vars, { model.getProps(state, 'sGmax') }];
                names = [names, {'sGmax'}];
                origin = [origin, class(model)];
                
            end
            
            if ~isempty(model.FacilityModel)
                
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
                
            end
        end

% ------------------------------------------------------------------------        
        
        function [eqs, state, src] = ...
                addBoundaryConditionsAndSources(model, eqs, names, types, ...
                                                state, forces)
            
            [p, s, mob, rho] = model.getProps(state, 'PhasePressures', ...
                                                     's'             , ...
                                                     'Mobility'     , ...
                                                     'Density');
            dissolved = {};
            comps = model.Components; 
                        
            [eqs, state, src] = addBoundaryConditionsAndSources@ReservoirModel(...
                model, eqs, names, types, state, p, s, mob, rho, dissolved, comps, forces);
        end

% ------------------------------------------------------------------------        
        
        function ctrl = validateDrivingForces(model, ctrl, index)
            ctrl = validateDrivingForces@ReservoirModel(model, ctrl, index);
        end

% ------------------------------------------------------------------------        
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ReservoirModel(model);
        end

% ------------------------------------------------------------------------        
        
        function state = initStateAD(model, state, vars, names, origin)

            toFacility = strcmp(origin, class(model.FacilityModel));
            remaining = ~toFacility;
            
            state = model.FacilityModel.initStateAD(state, vars(toFacility), ...
                                                    names(toFacility), ...
                                                    origin(toFacility));

            if model.hasInstantDissolution()
                                
                % The instant dissolution model has a special treatment of
                % saturation and dissolved CO2; both of which are stored in a
                % single primary variable 'X'.
                x_ix = strcmp(names, 'X');
                x = vars{x_ix};
                rs = model.getProp(state, 'rs');
                unsat = rs < model.fluid.dis_max;
                sg = x; sg(unsat) = 0;
                rs = x; rs(~unsat) = model.fluid.dis_max;
                state = model.setProp(state, 'saturation', {1-sg, sg});
                state = model.setProp(state, 'rs', rs);
                remaining(x_ix) = false;
            else
                % In models other than the instant dissolution model,
                % saturation is a primary variable
                sg_ix = strcmp(names, 'sG');
                sg = vars{sg_ix};
                state = model.setProp(state, 'saturation', {1-sg, sg});
                remaining(sg_ix) = false;
            end
            
            
            
            state = initStateAD@ReservoirModel(model, state, vars(remaining), ...
                                               names(remaining), ...
                                               origin(remaining));
        end

% ------------------------------------------------------------------------        
            
        function [v_eqs, tolerances, names] = ...
                getConvergenceValues(model, problem, varargin)
            [v_eqs, tolerances, names] = ...
                getConvergenceValues@ReservoirModel(model, problem, varargin{:});
        end

% ------------------------------------------------------------------------        
        
        function [model, state] = updateForChangedControls(model, state, forces)
            [model, state] = updateForChangedControls@ReservoirModel(model, state, forces);
        end

% ------------------------------------------------------------------------        

        function [eq, src] = ...
                addComponentContributions(model, cname, eq, component, src, force)     
            
            cnames = model.getComponentNames();
            ix = strcmpi(cnames, cname);

            % @@ Does not take dissolution into account.  Refine later.
            qC = src.phaseMass{find(ix)}; 
            
            if ~isempty(src.mapping)
                qC = src.mapping*qC;
            end
            
            cells = src.sourceCells;            
            eq(cells) = eq(cells) - qC;
            
        end
    
% ----------------------------------------------------------------------------
   
        function gdz = getGravityGradient(model)
            s  = model.operators;
            z = model.G.cells.z;
            if isfield(model.G, 'dz')
                z = z + model.G.dz;
            end
            gdz = model.gravity * s.Grad(z);
        end

% ----------------------------------------------------------------------------        

        function res = hasInstantDissolution(model)
            res = model.disgas && (model.fluid.dis_rate == 0 || ...
                                   model.fluid.dis_rate == Inf);
        end        
        
% ----------------------------------------------------------------------------        
        function res = hasRateDrivenDissolution(model)
            res = model.disgas && ~model.hasInstantDissolution();
        end        

% ----------------------------------------------------------------------------        
        function res = hasDissolution(model)
            res = hasInstantDissolution(model) || hasRateDrivenDissolution(model);
        end 
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
