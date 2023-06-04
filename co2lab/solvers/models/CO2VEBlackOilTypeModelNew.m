classdef CO2VEBlackOilTypeModelNew < ReservoirModel & GenericReservoirModel
    
    properties
        hysteresis
        disgas  % the name of this field should not be changed, in order to
                % play nice with the rest of AD-OO.
    end

    methods
        function model = CO2VEBlackOilTypeModelNew(Gt, rock2D, fluid, varargin)
           
            model = model@ReservoirModel(Gt, varargin{:});
            model.rock = rock2D;
            model.fluid = fluid;
            [model.water, model.oil, model.gas] = deal(true, false, true);
            model.gravity = norm(gravity);
            
            % determine if model includes hysteresis and/or dissolution
            model.hysteresis = model.fluid.res_gas > 0;
            model.disgas = isfield(model.fluid, 'dis_max'); 
            
            model = model.setupOperators(Gt, rock2D);
            model.Components = {MiscibleWaterComponent('brine', 1), ...
                                CO2Component('CO2', model.disgas, 2)};
        end
        
% ------------------------------------------------------------------------        

        function model = setupOperators(model, Gt, rock)
        % Compute vertially-integrated transmissibilities if not provided
            rock_tmp      = rock; 
            rock_tmp.perm = rock.perm .* Gt.cells.H; 
            T             = computeTrans(Gt, rock_tmp); 
            cf            = Gt.cells.faces(:, 1); 
            nf            = Gt.faces.num; 
            T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);       
            
            % Computing vertically-integrated pore - volume
            pv = poreVolume(Gt, rock); 
            pv = pv .* Gt.cells.H; 
            
            model.operators = setupOperatorsTPFA(Gt, rock, 'porv', pv, 'trans', T);
            model.operators.T_all = T;
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

            if model.hasRateDrivenDissolution()
                % Add conservation equation for dissolved gas
                [pv, pv0] = deal(model.getProp(state, 'PoreVolume'), ...
                                 model.getProp(state0, 'PoreVolume'));
                [rs, rs0] = deal(model.getProp(state, 'rs'), ...
                                 model.getProp(state0, 'rs'));
                [b, b0]   = deal(model.getProp(state, 'ShrinkageFactors'), ...
                                 model.getProp(state0, 'ShrinkageFactors'));
                [sW, sW0] = deal(model.getProp(state, 'sw'), ...
                                 model.getProp(state0, 'sw'));

                bW = b{1}; bW0 = b0{1}; 
                bG = b{2};
                
                acc_dis  = ((pv  .* bW  .* sW  .* rs) - ...
                            (pv0 .* bW0 .* sW0 .* rs0)) / dt;
                
                flux_dis = model.getProp(state, 'DissolvedFlux');
                
                eq_dis = model.operators.AccDiv(acc_dis, flux_dis);
                
                eta = model.computeDissolutionTransfer(state0, state, dt);
                eqs = [eqs, {eq_dis - eta}];
                names = [names, {'dissolution'}];
                types = [types, {'cell'}];
                
                % amount of saturation "eaten up" by dissolution.  Will  be
                % used to compute depletion of residual staturation further below.
                dSg = (eta * dt) ./ (pv .* bG); 
            else
                % @@ in the case of an instant dissolution model, setting dSg
                % to zero is not strictly correct.  Even if maximum
                % saturation is reached at the moment s > 0 in a cell, there
                % may still be later dissolution of CO2 from the fact that
                % brine flow may change the concentration of CO2 in the
                % column.  To account for this, we would need to do a strict
                % accounting of dissolved CO2 in each cell, just as we do in
                % the rate-driven case.
                dSg = 0;
            end
            
            
            if model.hysteresis
                sG = model.getProp(state, 'sg');
                sGmax = model.getProp(state, 'sGmax');
                sGmax0 = model.getProp(state0, 'sGmax');

                names = [names, {'hysteresis'}];
                types = [types, {'cell'}];
                
                fac = (1-model.fluid.res_water) / model.fluid.res_gas;
                eqs = [eqs, { sGmax - max(sGmax0 - dSg * fac, sG) }];
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
           
            flux_max_dis = model.getProp(state_rsmax, 'DissolvedFlux'); 

            % what is the maximum amount of CO2 that could be still absorbed by the
            % brine in the cell
            max_demand  = ((pv  .* bW  .* sW  .* rsmax) - ...
                           (pv0 .* bW0 .* sW0 .* rs0)) / dt;
            
            max_demand = model.operators.AccDiv(max_demand, flux_max_dis);
            
            % how much dissolved CO2 could maximally be supplied from the
            % remaining CO2 phase in the cell
            zeroADI = bW * 0; %@@
            max_supply = zeroADI + (pv0 .* bG0 .* sG0) / dt;
            
            %max_supply = max_supply - (pv .* bG .* sG) / dt;
            
            src = model.FacilityModel.getComponentSources(state);
            if ~isempty(src.cells)
                max_supply(src.cells) = max_supply(src.cells) + ...
                                        src.value{2} / model.fluid.rhoGS;
                % max_supply(src.cells) = max_supply(src.cells) + ...
                %                         max(src.value{2}, 0) / model.fluid.rhoGS;
            end
            
            phase_fluxes = model.getProp(state, 'ComponentPhaseFlux');
            co2_phase_flux = phase_fluxes{2, 2}; % second component in second phase
            co2_phase_flux = co2_phase_flux / model.fluid.rhoGS; % volume, not mass
            
            co2_inflow = -1 * model.operators.C' * co2_phase_flux;
            co2_inflow(co2_inflow < 0) = 0; % @@ necessary?
            max_supply = max_supply + value(co2_inflow);

            % The maximum amount of CO2 that can be dissolved in a cell
            % will be limited by the most strict of the two above bounds
            max_transfer = max(min(max_demand, max_supply), 0);
            
            % Newly drained areas with residual brine will reach max CO2
            % concentration instantly, regardless of rate.  We therefore set
            % a minimum rate
            rw = model.fluid.res_water;
            reswat_change = rw / (1 - rw) * (pv .* sG .* bW - pv0 .* sG0 .* bW0);
            
            min_transfer = zeroADI;
            sGmax = model.getProp(state, 'sGmax');
            imbibing = sG < sGmax;
            
            min_transfer(~imbibing) = ...
                model.fluid.dis_max * max(reswat_change(~imbibing), 0) / dt;
            
            % The amount of actually dissolved CO2 is also limited by the
            % actual dissolution rate
            eta = max(min(model.fluid.dis_rate, max_transfer), min_transfer);
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
            pvtprops = pvtprops.setStateFunction(...
                'ShrinkageFactors', ...
                ShrinkageFactors(model, 'usePhasePressures', false));
            
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
                [problem.primaryVariables, ix] = ...
                    model.stripVars(problem.primaryVariables, {'X'});
                dx(ix) = [];
            end
            
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                         drivingForces);
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
                bnd_flux_ixs = drivingForces.bc.face;
                bnd_fluxes = state.flux(bnd_flux_ixs, :);
                
                % add internal fluxes (@@ will overwrite boundary fluxes)
                qWG = model.getProp(state, 'PhaseFlux');
                state = storeFluxes(model, state, qWG{1}, [], qWG{2});
                
                % copy back the overwritten boundary fluxes
                state.flux(bnd_flux_ixs, :) = bnd_fluxes;
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

            qC = src.phaseMass{find(ix)}; % @@ Does not take dissolution into account.  Refine later.
            
            if ~isempty(src.mapping)
                qC = src.mapping*qC;
            end
            
            cells = src.sourceCells;            
            eq(cells) = eq(cells) - qC;
            
        end
    
% ----------------------------------------------------------------------------
   
        function gdz = getGravityGradient(model)
            s  = model.operators;
            Gt = model.G;
            gdz = model.gravity * s.Grad(Gt.cells.z);
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
        
    end
end

% function [problem, state] = getAdjointEquations(model, state0, state, dt, drivingForces, varargin)
