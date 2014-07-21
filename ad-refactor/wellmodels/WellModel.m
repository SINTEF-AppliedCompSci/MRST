classdef WellModel
    properties
        
        referencePressureIndex
        physicalModelIdentifier
        allowWellSignChange
        allowCrossflow
        allowControlSwitching
        verbose

        % Properties used during calculations
        bfactors
        surfaceDensities
        mobilities
        pressure
        referencePressure
        pseudocomponents
        maxPseudocomponents
        components
        nonlinearIteration
        W
    end
    
    methods
        function wmodel = WellModel()
            
        end
        
        function [sources, wellEqs, controlEqs, wellSol] = ...
                            computeWellFlux(wellmodel, model,...
                                            W, wellSol, bhp,...
                                            currentFluxes, pressure,...
                                            surfaceDensities, bfactors, ...
                                            mob, components, varargin)
                                        
            
            wellmodel.verbose = mrstVerbose();
            wellmodel.pseudocomponents = {{}};
            wellmodel.maxPseudocomponents = {{}};
            wellmodel.nonlinearIteration = nan;
            wellmodel.referencePressureIndex = 2;
            wellmodel.allowWellSignChange   = false;
            wellmodel.allowCrossflow        = true;
            wellmodel.allowControlSwitching = true;
            
            wellmodel = merge_options(wellmodel, varargin{:});
            
            % Store all current variables in the current well model
            % instance temporarily
            wellmodel.bfactors = bfactors;
            wellmodel.surfaceDensities = surfaceDensities;
            wellmodel.mobilities = mob;
            wellmodel.components = components;
            wellmodel.W = W;
            clear opt
            
            if isempty(W)
                sources = {};
                controlEqs = {};
                return
            end
            wc    = vertcat(W.cells);
            
            ncomp = numel(model.componentNames);
            ph = model.getActivePhases();
            nph = sum(ph);
            
            if ~iscell(pressure)
                % Support single reference pressure
                p = cell(ncomp, 1);
                [p{:}]  = deal(pressure);
            else
                p = pressure;
            end
            wellmodel.referencePressure = p{wellmodel.referencePressureIndex};
            wellmodel.pressure = p;
            clear pressure
            
            % Assertions to indirectly document the implementation
            assert(numel(p) == ncomp, ['Number of component pressure did not', ...
                                       ' match number of components!']);
            assert(numel(mob) == ncomp, ['Number of component mobilities did not', ...
                                       ' match number of components!']);
            assert(numel(currentFluxes) == nph, ['Number of phase fluxes did not', ...
                                                   ' match number of phases!']);
            assert(numel(components) == ncomp, ['Number of provided components did not', ...
                                                   ' match number of components!']);
            assert(numel(surfaceDensities) == nph && numel(bfactors) == nph, ...
                'Number of densities or bfactors did not match number of present phases!');
            
            % Update well pressure
            [wellSol, currentFluxes, bhp] = wellmodel.updatePressure(wellSol, currentFluxes, bhp, model);
            % Update well limits
            [wellSol, currentFluxes, bhp] = wellmodel.updateLimits(wellSol, currentFluxes, bhp, model);
            
            % Set up the actual equations
            [wellEqs, controlEqs, sources, wellSol] =...
                wellmodel.assembleEquations(wellSol, currentFluxes, bhp, model);
            
        end
        
        function [wellSol, q_s, bhp] = updateLimits(wellmodel, wellSol, q_s, bhp, model)
            if ~wellmodel.allowControlSwitching
                return
            end
            % First pad all values so that the well and fluxes are three
            % phase in appearance (with zeros for non-present values)
            [q_s, W_3ph, active] = padRatesAndCompi(q_s, wellmodel.W, model);
            
            [wellSol, withinLims] = updateSwitchedWellControls(wellmodel, model, wellSol, bhp, q_s);
            
            if all(withinLims)
                q_s = q_s(active);
                return
            end
            
            for k = 1:numel(wellSol)
                w = W_3ph(k);
                if ~withinLims(k)
                    v  = wellSol(k).val;
                    switch wellSol(k).type
                        case 'bhp'
                            bhp = assignValue(bhp,v,k);
                        case 'rate'
                            q_s{1} = assignValue(q_s{1}, v*w.compi(1), k);
                            q_s{2} = assignValue(q_s{2}, v*w.compi(2), k);
                            if numel(q_s)>2
                                q_s{3} = assignValue(q_s{3}, v*w.compi(1), k);
                            end
                        case 'orat'
                            q_s{2} = assignValue(q_s{2}, v , k);
                        case 'wrat'
                            q_s{1} = assignValue(q_s{1}, v , k);
                        case 'grat'
                            q_s{3} = assignValue(q_s{3}, v , k);
                    end % No good guess for qOs, etc...
                end
            end
            
            q_s = q_s(active);
        end
        
        function [wellSol, flux, bhp] = updatePressure(wellmodel, wellSol, flux, bhp, model)
            if isnan(wellmodel.nonlinearIteration)
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate welbore pressure-drop will never be updated']);
            end
            
            if wellmodel.nonlinearIteration == 1
                wellSol = updateConnectionDP(wellmodel, model, wellSol);
            end
            
            if all(isfinite(double(bhp)))
                return
            end
            % Otherwise we have some initialization to do
            pd   = double(wellmodel.referencePressure);
            for k = 1:numel(wellSol)
                if ~isfinite(wellSol(k).bhp)
                    v = pd(k) - wellSol(k).cdp(1);
                    if isempty(wellSol(k).sign)
                        dispif(wellmodel.verbose, 'No sign found in wellsol...\n');
                    else
                        v = v + 5*wellSol(k).sign*barsa;
                    end
                    bhp = assignValue(bhp, v, k);
                    wellSol(k).bhp = v;
                end
            end
        end
        
        
        
        function  [wellEqs, controlEqs, cq_s, wellSol] = assembleEquations(wellmodel,...
                                                wellSol, q_s, bhp, model)
            [wellEqs, cq_s, mix_s, status, cstatus] = ...
                            computeWellContributionsNew(wellmodel, model, wellSol, bhp, q_s);
            controlEqs =  setupWellControlEquations(wellSol, bhp, q_s, status, mix_s, model);
            
            
            % Update well properties which are not primary variables
            w = wellmodel.W;
            nConn       = cellfun(@numel, {w.cells})'; % # connections of each well
            perf2well   = rldecode((1:numel(w))', nConn);
            toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
            cq_sDb = cell2mat(toDouble(cq_s));
            
            for wnr = 1:numel(wellSol)
                ix = perf2well == wnr;
                wellSol(wnr).cqs     = cq_sDb(ix,:);
                wellSol(wnr).cstatus = cstatus(ix);
                wellSol(wnr).status  = status(wnr);
            end
        end
    end
    methods (Static)
        function id = getIdentifier(model)
            id = model.getPhaseNames();
            
            dg = isfield(model, 'disgas') && model.disgas;
            vo = isfield(model, 'vapoil') && model.vapoil;
            if dg || vo
                if (vo || dg) && isa(model, 'ThreePhaseBlackOilModel')
                    id = [id, 'VO'];
                elseif dg
                    id = [id, 'BO'];
                end
            end
        end
    end
end



