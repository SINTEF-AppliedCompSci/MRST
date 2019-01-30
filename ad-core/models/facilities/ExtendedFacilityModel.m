classdef ExtendedFacilityModel < FacilityModel
    properties
        
    end
    
    methods
        
        function n = getNumberOfComponents(fm)
            n = fm.ReservoirModel.getNumberOfComponents();
        end
        
        function n = getNumberOfPhases(fm)
            n = fm.ReservoirModel.getNumberOfPhases();
        end
        
        function src = getComponentSources(fm, state)
            [val, map] = fm.getProps(state, 'ComponentTotalFlux', 'FacilityWellMapping');
            src = struct('value', {val}, 'cells', map.cells);
        end
        
        function [eqs, names, types, state] = getModelEquations(facility, state0, state, dt, drivingForces)
            model = facility.ReservoirModel;
            [map, cflux] = facility.getProps(state, 'FacilityWellMapping', 'ComponentTotalFlux');
            [p, T] = facility.getSurfaceConditions();
            nph = model.getNumberOfPhases();
            surfaceRates = cell(1, nph);
            [surfaceRates{:}] = deal(0);
            
            wsum = map.perforationSum;
            for c = 1:numel(cflux)
                composition = model.Components{c}.getPhaseCompositionSurface(model, state, p, T);
                for ph = 1:nph
                    if ~isempty(composition{ph})
                        surfaceRates{ph} = surfaceRates{ph} + composition{ph}.*(wsum*cflux{c});
                    end
                end
            end
            rhoS = model.getSurfaceDensities();
            [eqs, names, types] = deal(cell(1, nph+1));
            
            % This is a temporary hack!
            q_s = state.FacilityState.primaryVariables(1:nph);
            bhp = state.FacilityState.primaryVariables{nph+1};
            [sn, phnames] = model.getPhaseNames();
            for ph = 1:nph
                surfaceRates{ph} = surfaceRates{ph}./rhoS(ph);
                eqs{ph} = surfaceRates{ph} - q_s{ph};
                names{ph} = [phnames{ph}, 'Wells'];
                types{ph} = 'perf';
            end
            
            mixs = value(surfaceRates);
            nact = numel(map.active);
            ctrl = cell(nact, 1);
            for i = 1:nact
                w = map.active(i);
                well = facility.WellModels{w};
                qs = cellfun(@(x) x(i), q_s, 'UniformOutput', false);
                ctrl{i} =  setupWellControlEquationsSingleWell(well, state0.wellSol(w), state.wellSol(w), bhp(i), qs, true, mixs(i, :), model);
            end
            eqs{end} = vertcat(ctrl{:});
            names{end} = 'closureWells';
            types{end} = 'well';
        end
        
        function state = applyWellLimits(fm, state)
            active = fm.getIndicesOfActiveWells(state.wellSol);
            for i = 1:numel(active)
                w = active(i);
                well = fm.WellModels{w};
                state.wellSol(w) = fm.applyWellLimitsWellSol(well, state.wellSol(w));
            end
        end

        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@FacilityModel(model, state, state0, dt, drivingForces);
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);
            
            if nw > 0
                rho = model.ReservoirModel.getProps(state, 'Density');
                rho = [rho{:}];
                [wc, perf2well] = model.getActiveWellCells(wellSol);
                rho = rho(wc, :);
                for i = 1:nw
                    wellNo = actWellIx(i);
                    wm = model.WellModels{wellNo};
                    % Possible index error for i here - should it be
                    % wellno?
                    rho_i = rho(perf2well == i, :);
                    wellSol(wellNo) = wm.updateConnectionPressureDropState(model.ReservoirModel, wellSol(wellNo), rho_i, rho_i);
                end
            end
            state.wellSol = wellSol;
            state = model.applyWellLimits(state);
        end
        
        function containers = getPropertyFunctions(model)
            containers = getPropertyFunctions@PhysicalModel(model);
            assert(not(isempty(model.FacilityFluxDiscretization)), ...
                'PropertyFunctions not initialized - did you call "validateModel"?');
            containers = [containers, {model.FacilityFluxDiscretization}];
        end
        
        function [p, T] = getSurfaceConditions(fm)
            p = 1*atm;
            T = 273.15 + 30;
        end
        
        function wellSol = applyWellLimitsWellSol(fm, well, wellSol)
            % Update solution variables and wellSol based on the well
            % limits. If limits have been reached, this function will
            % attempt to re-initialize the values and change the controls
            % so that the next step keeps within the prescribed ranges.
            withinLimits = true;
            if ~well.allowControlSwitching
                % We cannot change controls, so we return
                return
            end
            if isfield(well.W, 'status') && ~well.W.status
                % Well is inactive
                return
            end
            lims = well.W.lims;
            model = fm.ReservoirModel;

            if ~isnumeric(lims)
                phases = model.getPhaseNames();
            
                qs_t = 0;
                for i = 1:numel(phases)
                    qs_t = qs_t + wellSol.(['q', phases(i), 's']);
                end
                bhp = wellSol.bhp;
            
                if well.isInjector()
                    % Injectors have three possible limits:
                    % bhp:  Upper limit on pressure.
                    % rate: Upper limit on total surface rate.
                    % vrat: Lower limit on total surface rate.
                    modes   = {'bhp', 'rate', 'vrat'};
                    lims = well.setMissingLimits(lims, modes, inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is lower limit, switch default sign
                        lims.vrat = -inf;
                    end

                    flags = [value(bhp)  > lims.bhp, ...
                              qs_t       > lims.rate, ...
                              qs_t       < lims.vrat];
                else
                    % Producers have several possible limits:
                    % bhp:  Lower limit on pressure.
                    % orat: Lower limit on surface oil rate
                    % lrat: Lower limit on surface liquid (water + oil) rate
                    % grat: Lower limit on surface gas rate
                    % wrat: Lower limit on surface water rate
                    % vrat: Upper limit on total volumetric surface rate

                    modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                    lims = fm.setMissingLimits(lims, modes, -inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is upper limit, switch default sign
                        lims.vrat = inf;
                    end
                    [q_w, q_o, q_g, q_sl] = deal(0);
                    if model.water
                        q_w = wellSol.qWs;
                    end
                    if model.oil
                        q_o = wellSol.qOs;
                    end
                    if model.gas
                        q_g = wellSol.qGs;
                    end
                    if isprop(model, 'solvent') && model.solvent
                        q_sl = wellSol.qSs;
                    end
                    flags = [value(bhp) < lims.bhp,  ...
                        q_o          < lims.orat, ...
                        q_w + q_o    < lims.lrat, ...
                        q_g + q_sl   < lims.grat, ...
                        q_w          < lims.wrat, ...
                        qs_t         > lims.vrat];
                end
            else
                modes = {};
                flags = false;
                assert(isempty(lims) || isinf(lims))
            end
            % limits we need to check (all others than w.type):
            chkInx = ~strcmp(wellSol.type, modes);
            vltInx = find(flags(chkInx), 1);
            if ~isempty(vltInx)
                withinLimits = false;
                modes  = modes(chkInx);
                switchMode = modes{vltInx};
                fprintf('Well %s: Control mode changed from %s to %s.\n', wellSol.name, wellSol.type, switchMode);
                wellSol.type = switchMode;
                wellSol.val  = lims.(switchMode);
            end

            if ~withinLimits
                v  = wellSol.val;
                switch wellSol.type
                    case 'bhp'
                        wellSol.bhp = bhp;
                    case 'rate'
                        for ix = 1:numel(phases)
                            wellSol.(['q', phases(ix), 's']) = v*well.W.compi(ix);
                        end
                    case 'orat'
                        wellSol.qOs = v;
                    case 'wrat'
                        wellSol.qWs = v;
                    case 'grat'
                        wellSol.qGs = v;
                end % No good guess for qOs, etc...
            end
        end
        
        function lims = setMissingLimits(fm, lims, modes, val)
            missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
            for f = missing_fields
               lims = setfield(lims, f{:}, val);
            end
        end
    end
end