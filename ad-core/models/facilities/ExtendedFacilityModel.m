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
    end
end